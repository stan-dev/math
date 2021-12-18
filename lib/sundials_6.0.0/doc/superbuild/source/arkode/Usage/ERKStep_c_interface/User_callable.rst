.. ----------------------------------------------------------------
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

.. _ARKODE.Usage.ERKStep.UserCallable:

ERKStep User-callable functions
==================================

This section describes the functions that are called by the
user to setup and then solve an IVP using the ERKStep time-stepping
module. Some of these are required; however, starting with
:numref:`ARKODE.Usage.ERKStep.OptionalInputs`, the functions listed involve
optional inputs/outputs or restarting, and those paragraphs may be
skipped for a casual use of ARKODE's ERKStep module. In any case,
refer to the preceding section, :numref:`ARKODE.Usage.ERKStep.Skeleton`,
for the correct order of these calls.

On an error, each user-callable function returns a negative value  (or
``NULL`` if the function returns a pointer) and sends an error message
to the error handler routine, which prints the message to ``stderr``
by default. However, the user can set a file as error output or can
provide her own error handler function (see
:numref:`ARKODE.Usage.ERKStep.OptionalInputs` for details).



.. _ARKODE.Usage.ERKStep.Initialization:

ERKStep initialization and deallocation functions
------------------------------------------------------


.. c:function:: void* ERKStepCreate(ARKRhsFn f, realtype t0, N_Vector y0, SUNContext sunctx)

   This function allocates and initializes memory for a problem to
   be solved using the ERKStep time-stepping module in ARKODE.

   **Arguments:**
      * *f* -- the name of the C function (of type :c:func:`ARKRhsFn()`)
        defining the right-hand side function in
        :math:`\dot{y} = f(t,y)`.
      * *t0* -- the initial value of :math:`t`.
      * *y0* -- the initial condition vector :math:`y(t_0)`.
      * *sunctx* -- the :c:type:`SUNContext` object (see :numref:`SUNDIALS.SUNContext`)

   **Return value:**
      If successful, a pointer to initialized problem memory
      of type ``void*``, to be passed to all user-facing ERKStep routines
      listed below.  If unsuccessful, a ``NULL`` pointer will be
      returned, and an error message will be printed to ``stderr``.


.. c:function:: void ERKStepFree(void** arkode_mem)

   This function frees the problem memory *arkode_mem* created by
   :c:func:`ERKStepCreate`.

   **Arguments:**
      * *arkode_mem* -- pointer to the ERKStep memory block.

   **Return value:**  None



.. _ARKODE.Usage.ERKStep.Tolerances:

ERKStep tolerance specification functions
------------------------------------------------------

These functions specify the integration tolerances. One of them
**should** be called before the first call to
:c:func:`ERKStepEvolve()`; otherwise default values of ``reltol =
1e-4`` and ``abstol = 1e-9`` will be used, which may be entirely
incorrect for a specific problem.

The integration tolerances ``reltol`` and ``abstol`` define a vector
of error weights, ``ewt``.  In the case of
:c:func:`ERKStepSStolerances()`, this vector has components

.. code-block:: c

   ewt[i] = 1.0/(reltol*abs(y[i]) + abstol);

whereas in the case of :c:func:`ERKStepSVtolerances()` the vector components
are given by

.. code-block:: c

   ewt[i] = 1.0/(reltol*abs(y[i]) + abstol[i]);

This vector is used in all error tests, which use a weighted RMS norm
on all error-like vectors v:

.. math::
    \|v\|_{WRMS} = \left( \frac{1}{N} \sum_{i=1}^N (v_i\; ewt_i)^2 \right)^{1/2},

where :math:`N` is the problem dimension.

Alternatively, the user may supply a custom function to supply the
``ewt`` vector, through a call to :c:func:`ERKStepWFtolerances()`.



.. c:function:: int ERKStepSStolerances(void* arkode_mem, realtype reltol, realtype abstol)

   This function specifies scalar relative and absolute tolerances.

   **Arguments:**
      * *arkode_mem* -- pointer to the ERKStep memory block.
      * *reltol* -- scalar relative tolerance.
      * *abstol* -- scalar absolute tolerance.

   **Return value:**
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL*  if the ERKStep memory was ``NULL``
      * *ARK_NO_MALLOC*  if the ERKStep memory was not allocated by the time-stepping module
      * *ARK_ILL_INPUT* if an argument has an illegal value (e.g. a negative tolerance).



.. c:function:: int ERKStepSVtolerances(void* arkode_mem, realtype reltol, N_Vector abstol)

   This function specifies a scalar relative tolerance and a vector
   absolute tolerance (a potentially different absolute tolerance for
   each vector component).

   **Arguments:**
      * *arkode_mem* -- pointer to the ERKStep memory block.
      * *reltol* -- scalar relative tolerance.
      * *abstol* -- vector containing the absolute tolerances for each
        solution component.

   **Return value:**
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL*  if the ERKStep memory was ``NULL``
      * *ARK_NO_MALLOC*  if the ERKStep memory was not allocated by the time-stepping module
      * *ARK_ILL_INPUT* if an argument has an illegal value (e.g. a negative tolerance).



.. c:function:: int ERKStepWFtolerances(void* arkode_mem, ARKEwtFn efun)

   This function specifies a user-supplied function *efun* to compute
   the error weight vector ``ewt``.

   **Arguments:**
      * *arkode_mem* -- pointer to the ERKStep memory block.
      * *efun* -- the name of the function (of type :c:func:`ARKEwtFn()`)
        that implements the error weight vector computation.

   **Return value:**
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL*  if the ERKStep memory was ``NULL``
      * *ARK_NO_MALLOC*  if the ERKStep memory was not allocated by the time-stepping module




General advice on the choice of tolerances
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

For many users, the appropriate choices for tolerance values in
``reltol`` and ``abstol`` are a concern. The following pieces
of advice are relevant.

(1) The scalar relative tolerance ``reltol`` is to be set to control
    relative errors. So a value of :math:`10^{-4}` means that errors
    are controlled to .01%. We do not recommend using ``reltol`` larger
    than :math:`10^{-3}`. On the other hand, ``reltol`` should not be so
    small that it is comparable to the unit roundoff of the machine
    arithmetic (generally around :math:`10^{-15}` for double-precision).

(2) The absolute tolerances ``abstol`` (whether scalar or vector) need
    to be set to control absolute errors when any components of the
    solution vector :math:`y` may be so small that pure relative error
    control is meaningless.  For example, if :math:`y_i` starts at some
    nonzero value, but in time decays to zero, then pure relative
    error control on :math:`y_i` makes no sense (and is overly costly)
    after :math:`y_i` is below some noise level. Then ``abstol`` (if
    scalar) or ``abstol[i]`` (if a vector) needs to be set to that
    noise level. If the different components have different noise
    levels, then ``abstol`` should be a vector.  For example, see the
    example problem ``ark_robertson.c``, and the discussion
    of it in the ARKODE Examples Documentation :cite:p:`arkode_ex`.  In that
    problem, the three components vary between 0 and 1, and have
    different noise levels; hence the ``atols`` vector therein. It is
    impossible to give any general advice on ``abstol`` values,
    because the appropriate noise levels are completely
    problem-dependent. The user or modeler hopefully has some idea as
    to what those noise levels are.

(3) Finally, it is important to pick all the tolerance values
    conservatively, because they control the error committed on each
    individual step. The final (global) errors are an accumulation of
    those per-step errors, where that accumulation factor is
    problem-dependent.  A general rule of thumb is to reduce the
    tolerances by a factor of 10 from the actual desired limits on
    errors.  So if you want .01% relative accuracy (globally), a good
    choice for ``reltol`` is :math:`10^{-5}`.  In any case, it is
    a good idea to do a few experiments with the tolerances to see how
    the computed solution values vary as tolerances are reduced.



Advice on controlling nonphysical negative values
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In many applications, some components in the true solution are always
positive or non-negative, though at times very small.  In the
numerical solution, however, small negative (nonphysical) values
can then occur. In most cases, these values are harmless, and simply
need to be controlled, not eliminated, but in other cases any value
that violates a constraint may cause a simulation to halt. For both of
these scenarios the following pieces of advice are relevant.

(1) The best way to control the size of unwanted negative computed
    values is with tighter absolute tolerances.  Again this requires
    some knowledge of the noise level of these components, which may
    or may not be different for different components. Some
    experimentation may be needed.

(2) If output plots or tables are being generated, and it is important
    to avoid having negative numbers appear there (for the sake of
    avoiding a long explanation of them, if nothing else), then
    eliminate them, but only in the context of the output medium. Then
    the internal values carried by the solver are unaffected. Remember
    that a small negative value in :math:`y` returned by ERKStep, with
    magnitude comparable to ``abstol`` or less, is equivalent to zero
    as far as the computation is concerned.

(3) The user's right-hand side routine :math:`f`
    should never change a negative value in the solution vector :math:`y`
    to a non-negative value in attempt to "fix" this problem,
    since this can lead to numerical instability.  If the :math:`f`
    routine cannot tolerate a zero or negative value (e.g. because
    there is a square root or log), then the offending value should be
    changed to zero or a tiny positive number in a temporary variable
    (not in the input :math:`y` vector) for the purposes of computing
    :math:`f(t, y)`.

(4) ERKStep supports component-wise constraints on solution components,
    :math:`y_i < 0`, :math:`y_i \le 0`, , :math:`y_i > 0`, or
    :math:`y_i \ge 0`, through the user-callable function
    :c:func:`ERKStepSetConstraints`.  At each internal time step, if any
    constraint is violated then ERKStep will attempt a smaller time step
    that should not violate this constraint.  This reduced step size is
    chosen such that the step size is the largest possible but where the
    solution component satisfies the constraint.

(5) Positivity and non-negativity constraints on components can be
    enforced by use of the recoverable error return feature in the
    user-supplied right-hand side function, :math:`f`. When a
    recoverable error is encountered, ERKStep will retry the step with
    a smaller step size, which typically alleviates the problem.
    However, because this option involves some additional overhead
    cost, it should only be exercised if the use of absolute
    tolerances to control the computed values is unsuccessful.



.. _ARKODE.Usage.ERKStep.RootFinding:

Rootfinding initialization function
--------------------------------------

As described in :numref:`ARKODE.Mathematics.Rootfinding`, while
solving the IVP, ARKODE's time-stepping modules have the capability to
find the roots of a set of user-defined functions.  To activate the
root-finding algorithm, call the following function.  This is normally
called only once, prior to the first call to
:c:func:`ERKStepEvolve()`, but if the rootfinding problem is to be
changed during the solution, :c:func:`ERKStepRootInit()` can also be
called prior to a continuation call to :c:func:`ERKStepEvolve()`.


.. c:function:: int ERKStepRootInit(void* arkode_mem, int nrtfn, ARKRootFn g)

   Initializes a rootfinding problem to be solved during the
   integration of the ODE system.  It must be called after
   :c:func:`ERKStepCreate`, and before :c:func:`ERKStepEvolve()`.

   **Arguments:**
      * *arkode_mem* -- pointer to the ERKStep memory block.
      * *nrtfn* -- number of functions :math:`g_i`, an integer :math:`\ge` 0.
      * *g* -- name of user-supplied function, of type :c:func:`ARKRootFn()`,
        defining the functions :math:`g_i` whose roots are sought.

   **Return value:**
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL*  if the ERKStep memory was ``NULL``
      * *ARK_MEM_FAIL*  if there was a memory allocation failure
      * *ARK_ILL_INPUT* if *nrtfn* is greater than zero but *g* = ``NULL``.

   **Notes:**
      To disable the rootfinding feature after it has already
      been initialized, or to free memory associated with ERKStep's
      rootfinding module, call *ERKStepRootInit* with *nrtfn = 0*.

      Similarly, if a new IVP is to be solved with a call to
      :c:func:`ERKStepReInit()`, where the new IVP has no rootfinding
      problem but the prior one did, then call *ERKStepRootInit* with
      *nrtfn = 0*.




.. _ARKODE.Usage.ERKStep.Integration:

ERKStep solver function
-------------------------

This is the central step in the solution process -- the call to perform
the integration of the IVP.  One of the input arguments (*itask*)
specifies one of two modes as to where ERKStep is to return a
solution.  These modes are modified if the user has set a stop time
(with a call to the optional input function :c:func:`ERKStepSetStopTime()`) or
has requested rootfinding.



.. c:function:: int ERKStepEvolve(void* arkode_mem, realtype tout, N_Vector yout, realtype *tret, int itask)

   Integrates the ODE over an interval in :math:`t`.

   **Arguments:**
      * *arkode_mem* -- pointer to the ERKStep memory block.
      * *tout* -- the next time at which a computed solution is desired.
      * *yout* -- the computed solution vector.
      * *tret* -- the time corresponding to *yout* (output).
      * *itask* -- a flag indicating the job of the solver for the next
        user step.

        The *ARK_NORMAL* option causes the solver to take internal
        steps until it has just overtaken a user-specified output
        time, *tout*, in the direction of integration,
        i.e. :math:`t_{n-1} <` *tout* :math:`\le t_{n}` for forward
        integration, or :math:`t_{n} \le` *tout* :math:`< t_{n-1}` for
        backward integration.  It will then compute an approximation
        to the solution :math:`y(tout)` by interpolation (using one
        of the dense output routines described in
        :numref:`ARKODE.Mathematics.Interpolation`).

        The *ARK_ONE_STEP* option tells the solver to only take a
        single internal step :math:`y_{n-1} \to y_{n}` and then return
        control back to the calling program.  If this step will
        overtake *tout* then the solver will again return an
        interpolated result; otherwise it will return a copy of the
        internal solution :math:`y_{n}` in the vector *yout*.

   **Return value:**
      * *ARK_SUCCESS* if successful.
      * *ARK_ROOT_RETURN* if :c:func:`ERKStepEvolve()` succeeded, and
        found one or more roots.  If the number of root functions,
        *nrtfn*, is greater than 1, call
        :c:func:`ERKStepGetRootInfo()` to see which :math:`g_i` were
        found to have a root at (*\*tret*).
      * *ARK_TSTOP_RETURN* if :c:func:`ERKStepEvolve()` succeeded and
        returned at *tstop*.
      * *ARK_MEM_NULL* if the *arkode_mem* argument was ``NULL``.
      * *ARK_NO_MALLOC* if *arkode_mem* was not allocated.
      * *ARK_ILL_INPUT* if one of the inputs to
        :c:func:`ERKStepEvolve()` is illegal, or some other input to
        the solver was either illegal or missing.  Details will be
        provided in the error message.  Typical causes of this failure:

        (a) A component of the error weight vector became zero during
            internal time-stepping.

        (b) A root of one of the root functions was found both at a
            point :math:`t` and also very near :math:`t`.

        (c) The initial condition violates the inequality constraints.

      * *ARK_TOO_MUCH_WORK* if the solver took *mxstep* internal steps
        but could not reach *tout*.  The default value for *mxstep* is
        *MXSTEP_DEFAULT = 500*.
      * *ARK_TOO_MUCH_ACC* if the solver could not satisfy the accuracy
        demanded by the user for some internal step.
      * *ARK_ERR_FAILURE* if error test failures occurred either too many
        times (*ark_maxnef*) during one internal time step or occurred
        with :math:`|h| = h_{min}`.
      * *ARK_VECTOROP_ERR* a vector operation error occurred.

   **Notes:**
      The input vector *yout* can use the same memory as the
      vector *y0* of initial conditions that was passed to
      :c:func:`ERKStepCreate`.

      In *ARK_ONE_STEP* mode, *tout* is used only on the first call, and
      only to get the direction and a rough scale of the independent
      variable. All failure return values are negative and so testing the
      return argument for negative values will trap all
      :c:func:`ERKStepEvolve()` failures.

      Since interpolation may reduce the accuracy in the reported
      solution, if full method accuracy is desired the user should issue
      a call to :c:func:`ERKStepSetStopTime()` before the call to
      :c:func:`ERKStepEvolve()` to specify a fixed stop time to
      end the time step and return to the user.  Upon return from
      :c:func:`ERKStepEvolve()`, a copy of the internal solution
      :math:`y_{n}` will be returned in the vector *yout*.  Once the
      integrator returns at a *tstop* time, any future testing for
      *tstop* is disabled (and can be re-enabled only though a new call
      to :c:func:`ERKStepSetStopTime()`).

      On any error return in which one or more internal steps were taken
      by :c:func:`ERKStepEvolve()`, the returned values of *tret* and
      *yout* correspond to the farthest point reached in the integration.
      On all other error returns, *tret* and *yout* are left unchanged
      from those provided to the routine.




.. _ARKODE.Usage.ERKStep.OptionalInputs:

Optional input functions
-------------------------

There are numerous optional input parameters that control the behavior
of ERKStep, each of which may be modified from its default value through
calling an appropriate input function.  The following tables list all
optional input functions, grouped by which aspect of ERKStep they control.
Detailed information on the calling syntax and arguments for each
function are then provided following each table.

The optional inputs are grouped into the following categories:

* General ERKStep options (:numref:`ARKODE.Usage.ERKStep.ERKStepInputTable`),

* IVP method solver options (:numref:`ARKODE.Usage.ERKStep.ERKStepMethodInputTable`),

* Step adaptivity solver options (:numref:`ARKODE.Usage.ERKStep.ERKStepAdaptivityInputTable`), and

* Rootfinding options (:numref:`ARKODE.Usage.ERKStep.ERKStepRootfindingInputTable`).

For the most casual use of ERKStep, relying on the default set of
solver parameters, the reader can skip to section on user-supplied
functions, :numref:`ARKODE.Usage.UserSupplied`.

We note that, on an error return, all of the optional input functions send an
error message to the error handler function. All error return values are
negative, so a test on the return arguments for negative values will catch all
errors. Finally, a call to an ``ERKStepSet***`` function can generally be made
from the user's calling program at any time and, if successful, takes effect
immediately. ``ERKStepSet***`` functions that cannot be called at any time note
this in the "**Notes**:" section of the function documentation.



.. _ARKODE.Usage.ERKStep.ERKStepInput:

Optional inputs for ERKStep
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. _ARKODE.Usage.ERKStep.ERKStepInputTable:
.. table:: Optional inputs for ERKStep

   +----------------------------------------------------+-----------------------------------------+------------------------+
   | Optional input                                     | Function name                           |  Default               |
   +----------------------------------------------------+-----------------------------------------+------------------------+
   | Return ERKStep solver parameters to their defaults | :c:func:`ERKStepSetDefaults()`          |  internal              |
   +----------------------------------------------------+-----------------------------------------+------------------------+
   | Set dense output interpolation type                | :c:func:`ERKStepSetInterpolantType()`   | ``ARK_INTERP_HERMITE`` |
   +----------------------------------------------------+-----------------------------------------+------------------------+
   | Set dense output polynomial degree                 | :c:func:`ERKStepSetInterpolantDegree()` |  5                     |
   +----------------------------------------------------+-----------------------------------------+------------------------+
   | Supply a pointer to a diagnostics output file      | :c:func:`ERKStepSetDiagnostics()`       | ``NULL``               |
   +----------------------------------------------------+-----------------------------------------+------------------------+
   | Supply a pointer to an error output file           | :c:func:`ERKStepSetErrFile()`           | ``stderr``             |
   +----------------------------------------------------+-----------------------------------------+------------------------+
   | Supply a custom error handler function             | :c:func:`ERKStepSetErrHandlerFn()`      |  internal fn           |
   +----------------------------------------------------+-----------------------------------------+------------------------+
   | Disable time step adaptivity (fixed-step mode)     | :c:func:`ERKStepSetFixedStep()`         |  disabled              |
   +----------------------------------------------------+-----------------------------------------+------------------------+
   | Supply an initial step size to attempt             | :c:func:`ERKStepSetInitStep()`          |  estimated             |
   +----------------------------------------------------+-----------------------------------------+------------------------+
   | Maximum no. of warnings for :math:`t_n+h = t_n`    | :c:func:`ERKStepSetMaxHnilWarns()`      |  10                    |
   +----------------------------------------------------+-----------------------------------------+------------------------+
   | Maximum no. of internal steps before *tout*        | :c:func:`ERKStepSetMaxNumSteps()`       |  500                   |
   +----------------------------------------------------+-----------------------------------------+------------------------+
   | Maximum absolute step size                         | :c:func:`ERKStepSetMaxStep()`           | :math:`\infty`         |
   +----------------------------------------------------+-----------------------------------------+------------------------+
   | Minimum absolute step size                         | :c:func:`ERKStepSetMinStep()`           |  0.0                   |
   +----------------------------------------------------+-----------------------------------------+------------------------+
   | Set a value for :math:`t_{stop}`                   | :c:func:`ERKStepSetStopTime()`          | :math:`\infty`         |
   +----------------------------------------------------+-----------------------------------------+------------------------+
   | Supply a pointer for user data                     | :c:func:`ERKStepSetUserData()`          | ``NULL``               |
   +----------------------------------------------------+-----------------------------------------+------------------------+
   | Maximum no. of ERKStep error test failures         | :c:func:`ERKStepSetMaxErrTestFails()`   |  7                     |
   +----------------------------------------------------+-----------------------------------------+------------------------+
   | Set inequality constraints on solution             | :c:func:`ERKStepSetConstraints()`       | ``NULL``               |
   +----------------------------------------------------+-----------------------------------------+------------------------+
   | Set max number of constraint failures              | :c:func:`ERKStepSetMaxNumConstrFails()` |  10                    |
   +----------------------------------------------------+-----------------------------------------+------------------------+



.. c:function:: int ERKStepSetDefaults(void* arkode_mem)

   Resets all optional input parameters to ERKStep's original
   default values.

   **Arguments:**
      * *arkode_mem* -- pointer to the ERKStep memory block.

   **Return value:**
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL* if the ERKStep memory is ``NULL``
      * *ARK_ILL_INPUT* if an argument has an illegal value

   **Notes:**
      Does not change problem-defining function pointer *f*
      or the *user_data* pointer.

      Also leaves alone any data structures or options related to
      root-finding (those can be reset using :c:func:`ERKStepRootInit()`).



.. c:function:: int ERKStepSetInterpolantType(void* arkode_mem, int itype)

   Specifies use of the Lagrange or Hermite interpolation modules (used for
   dense output -- interpolation of solution output values and implicit
   method predictors).

   **Arguments:**
      * *arkode_mem* -- pointer to the ERKStep memory block.
      * *itype* -- requested interpolant type (``ARK_INTERP_HERMITE`` or ``ARK_INTERP_LAGRANGE``)

   **Return value:**
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL* if the ERKStep memory is ``NULL``
      * *ARK_MEM_FAIL* if the interpolation module cannot be allocated
      * *ARK_ILL_INPUT* if the *itype* argument is not recognized or the
        interpolation module has already been initialized

   **Notes:**
      The Hermite interpolation module is described in
      :numref:`ARKODE.Mathematics.Interpolation.Hermite`, and the Lagrange interpolation module
      is described in :numref:`ARKODE.Mathematics.Interpolation.Lagrange`.

      This routine frees any previously-allocated interpolation module, and re-creates
      one according to the specified argument.  Thus any previous calls to
      :c:func:`ERKStepSetInterpolantDegree()` will be nullified.

      This routine must be called *after* the call to :c:func:`ERKStepCreate`.
      After the first call to :c:func:`ERKStepEvolve()` the interpolation type may
      not be changed without first calling :c:func:`ERKStepReInit()`.

      If this routine is not called, the Hermite interpolation module will be used.



.. c:function:: int ERKStepSetInterpolantDegree(void* arkode_mem, int degree)

   Specifies the degree of the polynomial interpolant
   used for dense output (i.e. interpolation of solution output values
   and implicit method predictors).

   **Arguments:**
      * *arkode_mem* -- pointer to the ERKStep memory block.
      * *degree* -- requested polynomial degree.

   **Return value:**
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL* if the ERKStep memory or interpolation module are ``NULL``
      * *ARK_INTERP_FAIL* if this is called after :c:func:`ERKStepEvolve()`
      * *ARK_ILL_INPUT* if an argument has an illegal value or the
        interpolation module has already been initialized

   **Notes:**
      Allowed values are between 0 and 5.

      This routine should be called *after* :c:func:`ERKStepCreate` and *before*
      :c:func:`ERKStepEvolve()`. After the first call to :c:func:`ERKStepEvolve()`
      the interpolation degree may not be changed without first calling
      :c:func:`ERKStepReInit()`.

      If a user calls both this routine and :c:func:`ERKStepSetInterpolantType()`, then
      :c:func:`ERKStepSetInterpolantType()` must be called first.

      Since the accuracy of any polynomial interpolant is limited by the accuracy of
      the time-step solutions on which it is based, the *actual* polynomial degree that
      is used by ERKStep will be the minimum of :math:`q-1` and the input *degree*,
      where :math:`q` is the order of accuracy for the time integration method.



.. c:function:: int ERKStepSetDenseOrder(void* arkode_mem, int dord)

   *This function is deprecated, and will be removed in a future release.
   Users should transition to calling* :c:func:`ERKStepSetInterpolantDegree()`
   *instead.*



.. c:function:: int ERKStepSetDiagnostics(void* arkode_mem, FILE* diagfp)

   Specifies the file pointer for a diagnostics file where
   all ERKStep step adaptivity and solver information is written.

   **Arguments:**
      * *arkode_mem* -- pointer to the ERKStep memory block.
      * *diagfp* -- pointer to the diagnostics output file.

   **Return value:**
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL* if the ERKStep memory is ``NULL``
      * *ARK_ILL_INPUT* if an argument has an illegal value

   **Notes:**
      This parameter can be ``stdout`` or ``stderr``, although the
      suggested approach is to specify a pointer to a unique file opened
      by the user and returned by ``fopen``.  If not called, or if called
      with a ``NULL`` file pointer, all diagnostics output is disabled.

      When run in parallel, only one process should set a non-NULL value
      for this pointer, since statistics from all processes would be
      identical.



.. c:function:: int ERKStepSetErrFile(void* arkode_mem, FILE* errfp)

   Specifies a pointer to the file where all ERKStep warning and error
   messages will be written if the default internal error handling
   function is used.

   **Arguments:**
      * *arkode_mem* -- pointer to the ERKStep memory block.
      * *errfp* -- pointer to the output file.

   **Return value:**
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL* if the ERKStep memory is ``NULL``
      * *ARK_ILL_INPUT* if an argument has an illegal value

   **Notes:**
      The default value for *errfp* is ``stderr``.

      Passing a ``NULL`` value disables all future error message output
      (except for the case wherein the ERKStep memory pointer is
      ``NULL``).  This use of the function is strongly discouraged.

      If used, this routine should be called before any other
      optional input functions, in order to take effect for subsequent
      error messages.



.. c:function:: int ERKStepSetErrHandlerFn(void* arkode_mem, ARKErrHandlerFn ehfun, void* eh_data)

   Specifies the optional user-defined function to be used
   in handling error messages.

   **Arguments:**
      * *arkode_mem* -- pointer to the ERKStep memory block.
      * *ehfun* -- name of user-supplied error handler function.
      * *eh_data* -- pointer to user data passed to *ehfun* every time
        it is called.

   **Return value:**
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL* if the ERKStep memory is ``NULL``
      * *ARK_ILL_INPUT* if an argument has an illegal value

   **Notes:**
      Error messages indicating that the ERKStep solver memory is
      ``NULL`` will always be directed to ``stderr``.




.. c:function:: int ERKStepSetFixedStep(void* arkode_mem, realtype hfixed)

   Disabled time step adaptivity within ERKStep, and specifies the
   fixed time step size to use for the following internal step(s).

   **Arguments:**
      * *arkode_mem* -- pointer to the ERKStep memory block.
      * *hfixed* -- value of the fixed step size to use.

   **Return value:**
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL* if the ERKStep memory is ``NULL``
      * *ARK_ILL_INPUT* if an argument has an illegal value

   **Notes:**
      Pass 0.0 to return ERKStep to the default (adaptive-step) mode.

      Use of this function is not generally recommended, since we it gives no
      assurance of the validity of the computed solutions.  It is
      primarily provided for code-to-code verification testing purposes.

      When using :c:func:`ERKStepSetFixedStep()`, any values provided to
      the functions
      :c:func:`ERKStepSetInitStep()`,
      :c:func:`ERKStepSetAdaptivityFn()`,
      :c:func:`ERKStepSetMaxErrTestFails()`,
      :c:func:`ERKStepSetAdaptivityMethod()`,
      :c:func:`ERKStepSetCFLFraction()`,
      :c:func:`ERKStepSetErrorBias()`,
      :c:func:`ERKStepSetFixedStepBounds()`,
      :c:func:`ERKStepSetMaxEFailGrowth()`,
      :c:func:`ERKStepSetMaxFirstGrowth()`,
      :c:func:`ERKStepSetMaxGrowth()`,
      :c:func:`ERKStepSetMinReduction()`,
      :c:func:`ERKStepSetSafetyFactor()`,
      :c:func:`ERKStepSetSmallNumEFails()` and
      :c:func:`ERKStepSetStabilityFn()`
      will be ignored, since temporal adaptivity is disabled.

      If both :c:func:`ERKStepSetFixedStep()` and
      :c:func:`ERKStepSetStopTime()` are used, then the fixed step size
      will be used for all steps until the final step preceding the
      provided stop time (which may be shorter).  To resume use of the
      previous fixed step size, another call to
      :c:func:`ERKStepSetFixedStep()` must be made prior to calling
      :c:func:`ERKStepEvolve()` to resume integration.

      It is *not* recommended that :c:func:`ERKStepSetFixedStep()` be used
      in concert with :c:func:`ERKStepSetMaxStep()` or
      :c:func:`ERKStepSetMinStep()`, since at best those latter two
      routines will provide no useful information to the solver, and at
      worst they may interfere with the desired fixed step size.




.. c:function:: int ERKStepSetInitStep(void* arkode_mem, realtype hin)

   Specifies the initial time step size ERKStep should use after
   initialization, re-initialization, or resetting.

   **Arguments:**
      * *arkode_mem* -- pointer to the ERKStep memory block.
      * *hin* -- value of the initial step to be attempted :math:`(\ne 0)`.

   **Return value:**
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL* if the ERKStep memory is ``NULL``
      * *ARK_ILL_INPUT* if an argument has an illegal value

   **Notes:**
      Pass 0.0 to use the default value.

      By default, ERKStep estimates the initial step size to be
      :math:`h = \sqrt{\dfrac{2}{\left\| \ddot{y} \right\|}}`, where
      :math:`\ddot{y}` is an estimate of the second derivative of the
      solution at :math:`t_0`.

      This routine will also reset the step size and error history.



.. c:function:: int ERKStepSetMaxHnilWarns(void* arkode_mem, int mxhnil)

   Specifies the maximum number of messages issued by the
   solver to warn that :math:`t+h=t` on the next internal step, before
   ERKStep will instead return with an error.

   **Arguments:**
      * *arkode_mem* -- pointer to the ERKStep memory block.
      * *mxhnil* -- maximum allowed number of warning messages :math:`(>0)`.

   **Return value:**
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL* if the ERKStep memory is ``NULL``
      * *ARK_ILL_INPUT* if an argument has an illegal value

   **Notes:**
      The default value is 10; set *mxhnil* to zero to specify
      this default.

      A negative value indicates that no warning messages should be issued.




.. c:function:: int ERKStepSetMaxNumSteps(void* arkode_mem, long int mxsteps)

   Specifies the maximum number of steps to be taken by the
   solver in its attempt to reach the next output time, before ERKStep
   will return with an error.

   **Arguments:**
      * *arkode_mem* -- pointer to the ERKStep memory block.
      * *mxsteps* -- maximum allowed number of internal steps.

   **Return value:**
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL* if the ERKStep memory is ``NULL``
      * *ARK_ILL_INPUT* if an argument has an illegal value

   **Notes:**
      Passing *mxsteps* = 0 results in ERKStep using the
      default value (500).

      Passing *mxsteps* < 0 disables the test (not recommended).



.. c:function:: int ERKStepSetMaxStep(void* arkode_mem, realtype hmax)

   Specifies the upper bound on the magnitude of the time step size.

   **Arguments:**
      * *arkode_mem* -- pointer to the ERKStep memory block.
      * *hmax* -- maximum absolute value of the time step size :math:`(\ge 0)`.

   **Return value:**
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL* if the ERKStep memory is ``NULL``
      * *ARK_ILL_INPUT* if an argument has an illegal value

   **Notes:**
      Pass *hmax* :math:`\le 0.0` to set the default value of :math:`\infty`.



.. c:function:: int ERKStepSetMinStep(void* arkode_mem, realtype hmin)

   Specifies the lower bound on the magnitude of the time step size.

   **Arguments:**
      * *arkode_mem* -- pointer to the ERKStep memory block.
      * *hmin* -- minimum absolute value of the time step size :math:`(\ge 0)`.

   **Return value:**
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL* if the ERKStep memory is ``NULL``
      * *ARK_ILL_INPUT* if an argument has an illegal value

   **Notes:**
      Pass *hmin* :math:`\le 0.0` to set the default value of 0.



.. c:function:: int ERKStepSetStopTime(void* arkode_mem, realtype tstop)

   Specifies the value of the independent variable
   :math:`t` past which the solution is not to proceed.

   **Arguments:**
      * *arkode_mem* -- pointer to the ERKStep memory block.
      * *tstop* -- stopping time for the integrator.

   **Return value:**
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL* if the ERKStep memory is ``NULL``
      * *ARK_ILL_INPUT* if an argument has an illegal value

   **Notes:**
      The default is that no stop time is imposed.




.. c:function:: int ERKStepSetUserData(void* arkode_mem, void* user_data)

   Specifies the user data block *user_data* and
   attaches it to the main ERKStep memory block.

   **Arguments:**
      * *arkode_mem* -- pointer to the ERKStep memory block.
      * *user_data* -- pointer to the user data.

   **Return value:**
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL* if the ERKStep memory is ``NULL``
      * *ARK_ILL_INPUT* if an argument has an illegal value

   **Notes:**
      If specified, the pointer to *user_data* is passed to all
      user-supplied functions for which it is an argument; otherwise
      ``NULL`` is passed.




.. c:function:: int ERKStepSetMaxErrTestFails(void* arkode_mem, int maxnef)

   Specifies the maximum number of error test failures
   permitted in attempting one step, before returning with an error.

   **Arguments:**
      * *arkode_mem* -- pointer to the ERKStep memory block.
      * *maxnef* -- maximum allowed number of error test failures :math:`(>0)`.

   **Return value:**
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL* if the ERKStep memory is ``NULL``
      * *ARK_ILL_INPUT* if an argument has an illegal value

   **Notes:**
      The default value is 7; set *maxnef* :math:`\le 0`
      to specify this default.



.. c:function:: int ERKStepSetConstraints(void* arkode_mem, N_Vector constraints)

   Specifies a vector defining inequality constraints for each component of the
   solution vector :math:`y`.

   **Arguments:**
      * *arkode_mem* -- pointer to the ERKStep memory block.
      * *constraints* -- vector of constraint flags. Each component specifies
        the type of solution constraint:

        .. math::

           \texttt{constraints[i]} = \left\{ \begin{array}{rcl}
           0.0  &\Rightarrow\;& \text{no constraint is imposed on}\; y_i,\\
           1.0  &\Rightarrow\;& y_i \geq 0,\\
           -1.0  &\Rightarrow\;& y_i \leq 0,\\
           2.0  &\Rightarrow\;& y_i > 0,\\
           -2.0  &\Rightarrow\;& y_i < 0.\\
           \end{array}\right.

   **Return value:**
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL* if the ERKStep memory is ``NULL``
      * *ARK_ILL_INPUT* if the constraints vector contains illegal values

   **Notes:**
      The presence of a non-``NULL`` constraints vector that is not 0.0
      in all components will cause constraint checking to be performed. However, a
      call with 0.0 in all components of ``constraints`` will result in an illegal
      input return. A ``NULL`` constraints vector will disable constraint checking.

      After a call to :c:func:`ERKStepResize()` inequality constraint checking
      will be disabled and a call to :c:func:`ERKStepSetConstraints()` is
      required to re-enable constraint checking.

      Since constraint-handling is performed through cutting time steps that would
      violate the constraints, it is possible that this feature will cause some
      problems to fail due to an inability to enforce constraints even at the
      minimum time step size.  Additionally, the features :c:func:`ERKStepSetConstraints()`
      and :c:func:`ERKStepSetFixedStep()` are incompatible, and should not be used
      simultaneously.


.. c:function:: int ERKStepSetMaxNumConstrFails(void* arkode_mem, int maxfails)

   Specifies the maximum number of constraint failures in a step before ERKStep
   will return with an error.

   **Arguments:**
      * *arkode_mem* -- pointer to the ERKStep memory block.
      * *maxfails* -- maximum allowed number of constrain failures.

   **Return value:**
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL* if the ERKStep memory is ``NULL``

   **Notes:**
      Passing *maxfails* <= 0 results in ERKStep using the
      default value (10).



.. _ARKODE.Usage.ERKStep.ERKStepMethodInput:

Optional inputs for IVP method selection
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. _ARKODE.Usage.ERKStep.ERKStepMethodInputTable:
.. table:: Optional inputs for IVP method selection

   +----------------------------------+---------------------------------+------------------+
   | Optional input                   | Function name                   | Default          |
   +----------------------------------+---------------------------------+------------------+
   | Set integrator method order      | :c:func:`ERKStepSetOrder()`     | 4                |
   +----------------------------------+---------------------------------+------------------+
   | Set explicit RK table            | :c:func:`ERKStepSetTable()`     | internal         |
   +----------------------------------+---------------------------------+------------------+
   | Specify explicit RK table number | :c:func:`ERKStepSetTableNum()`  | internal         |
   +----------------------------------+---------------------------------+------------------+



.. c:function:: int ERKStepSetOrder(void* arkode_mem, int ord)

   Specifies the order of accuracy for the ERK integration method.

   **Arguments:**
      * *arkode_mem* -- pointer to the ERKStep memory block.
      * *ord* -- requested order of accuracy.

   **Return value:**
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL* if the ERKStep memory is ``NULL``
      * *ARK_ILL_INPUT* if an argument has an illegal value

   **Notes:**
      The allowed values are :math:`2 \le` *ord* :math:`\le
      8`.  Any illegal input will result in the default value of 4.

      Since *ord* affects the memory requirements for the internal
      ERKStep memory block, it cannot be changed after the first call to
      :c:func:`ERKStepEvolve()`, unless :c:func:`ERKStepReInit()` is called.



.. c:function:: int ERKStepSetTable(void* arkode_mem, ARKodeButcherTable B)

   Specifies a customized Butcher table for the ERK method.

   **Arguments:**
      * *arkode_mem* -- pointer to the ERKStep memory block.
      * *B* -- the Butcher table for the explicit RK method.

   **Return value:**
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL* if the ERKStep memory is ``NULL``
      * *ARK_ILL_INPUT* if an argument has an illegal value

   **Notes:**

      For a description of the :c:type:`ARKodeButcherTable` type and related
      functions for creating Butcher tables, see :numref:`ARKodeButcherTable`.

      No error checking is performed to ensure that either the method order *p* or
      the embedding order *q* specified in the Butcher table structure correctly
      describe the coefficients in the Butcher table.

      Error checking is performed to ensure that the Butcher table is strictly
      lower-triangular (i.e. that it specifies an ERK method).

      If the Butcher table does not contain an embedding, the user *must* call
      :c:func:`ERKStepSetFixedStep()` to enable fixed-step mode and set the desired
      time step size.



.. c:function:: int ERKStepSetTableNum(void* arkode_mem, ARKODE_ERKTableID etable)

   Indicates to use a specific built-in Butcher table for the ERK method.

   **Arguments:**
      * *arkode_mem* -- pointer to the ERKStep memory block.
      * *etable* -- index of the Butcher table.

   **Return value:**
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL* if the ERKStep memory is ``NULL``
      * *ARK_ILL_INPUT* if an argument has an illegal value

   **Notes:**
      *etable* should match an existing explicit method from
      :numref:`Butcher.explicit`.  Error-checking is performed
      to ensure that the table exists, and is not implicit.







.. _ARKODE.Usage.ERKStep.ERKStepAdaptivityInput:

Optional inputs for time step adaptivity
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The mathematical explanation of ARKODE's time step adaptivity
algorithm, including how each of the parameters below is used within
the code, is provided in :numref:`ARKODE.Mathematics.Adaptivity`.


.. _ARKODE.Usage.ERKStep.ERKStepAdaptivityInputTable:
.. table:: Optional inputs for time step adaptivity

   +-----------------------------------------------------------+----------------------------------------+-----------+
   | Optional input                                            | Function name                          | Default   |
   +-----------------------------------------------------------+----------------------------------------+-----------+
   | Set a custom time step adaptivity function                | :c:func:`ERKStepSetAdaptivityFn()`     | internal  |
   +-----------------------------------------------------------+----------------------------------------+-----------+
   | Choose an existing time step adaptivity method            | :c:func:`ERKStepSetAdaptivityMethod()` | 0         |
   +-----------------------------------------------------------+----------------------------------------+-----------+
   | Explicit stability safety factor                          | :c:func:`ERKStepSetCFLFraction()`      | 0.5       |
   +-----------------------------------------------------------+----------------------------------------+-----------+
   | Time step error bias factor                               | :c:func:`ERKStepSetErrorBias()`        | 1.5       |
   +-----------------------------------------------------------+----------------------------------------+-----------+
   | Bounds determining no change in step size                 | :c:func:`ERKStepSetFixedStepBounds()`  | 1.0  1.5  |
   +-----------------------------------------------------------+----------------------------------------+-----------+
   | Maximum step growth factor on error test fail             | :c:func:`ERKStepSetMaxEFailGrowth()`   | 0.3       |
   +-----------------------------------------------------------+----------------------------------------+-----------+
   | Maximum first step growth factor                          | :c:func:`ERKStepSetMaxFirstGrowth()`   | 10000.0   |
   +-----------------------------------------------------------+----------------------------------------+-----------+
   | Maximum allowed general step growth factor                | :c:func:`ERKStepSetMaxGrowth()`        | 20.0      |
   +-----------------------------------------------------------+----------------------------------------+-----------+
   | Minimum allowed step reduction factor on error test fail  | :c:func:`ERKStepSetMinReduction()`     | 0.1       |
   +-----------------------------------------------------------+----------------------------------------+-----------+
   | Time step safety factor                                   | :c:func:`ERKStepSetSafetyFactor()`     | 0.96      |
   +-----------------------------------------------------------+----------------------------------------+-----------+
   | Error fails before MaxEFailGrowth takes effect            | :c:func:`ERKStepSetSmallNumEFails()`   | 2         |
   +-----------------------------------------------------------+----------------------------------------+-----------+
   | Explicit stability function                               | :c:func:`ERKStepSetStabilityFn()`      | none      |
   +-----------------------------------------------------------+----------------------------------------+-----------+



.. c:function:: int ERKStepSetAdaptivityFn(void* arkode_mem, ARKAdaptFn hfun, void* h_data)

   Sets a user-supplied time-step adaptivity function.

   **Arguments:**
      * *arkode_mem* -- pointer to the ERKStep memory block.
      * *hfun* -- name of user-supplied adaptivity function.
      * *h_data* -- pointer to user data passed to *hfun* every time
        it is called.

   **Return value:**
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL* if the ERKStep memory is ``NULL``
      * *ARK_ILL_INPUT* if an argument has an illegal value

   **Notes:**
      This function should focus on accuracy-based time step
      estimation; for stability based time steps the function
      :c:func:`ERKStepSetStabilityFn()` should be used instead.



.. c:function:: int ERKStepSetAdaptivityMethod(void* arkode_mem, int imethod, int idefault, int pq, realtype* adapt_params)

   Specifies the method (and associated parameters) used for time step adaptivity.

   **Arguments:**
      * *arkode_mem* -- pointer to the ERKStep memory block.
      * *imethod* -- accuracy-based adaptivity method choice
        (0 :math:`\le` `imethod` :math:`\le` 5):
        0 is PID, 1 is PI, 2 is I, 3 is explicit Gustafsson, 4 is
        implicit Gustafsson, and 5 is the ImEx Gustafsson.
      * *idefault* -- flag denoting whether to use default adaptivity
        parameters (1), or that they will be supplied in the
        *adapt_params* argument (0).
      * *pq* -- flag denoting whether to use the embedding order of
        accuracy *p* (0) or the method order of accuracy *q* (1)
        within the adaptivity algorithm.  *p* is the default.
      * *adapt_params[0]* -- :math:`k_1` parameter within accuracy-based adaptivity algorithms.
      * *adapt_params[1]* -- :math:`k_2` parameter within accuracy-based adaptivity algorithms.
      * *adapt_params[2]* -- :math:`k_3` parameter within accuracy-based adaptivity algorithms.

   **Return value:**
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL* if the ERKStep memory is ``NULL``
      * *ARK_ILL_INPUT* if an argument has an illegal value

   **Notes:**
      If custom parameters are supplied, they will be checked
      for validity against published stability intervals.  If other
      parameter values are desired, it is recommended to instead provide
      a custom function through a call to :c:func:`ERKStepSetAdaptivityFn()`.



.. c:function:: int ERKStepSetCFLFraction(void* arkode_mem, realtype cfl_frac)

   Specifies the fraction of the estimated explicitly stable step to use.

   **Arguments:**
      * *arkode_mem* -- pointer to the ERKStep memory block.
      * *cfl_frac* -- maximum allowed fraction of explicitly stable step (default is 0.5).

   **Return value:**
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL* if the ERKStep memory is ``NULL``
      * *ARK_ILL_INPUT* if an argument has an illegal value

   **Notes:**
      Any non-positive parameter will imply a reset to the default
      value.



.. c:function:: int ERKStepSetErrorBias(void* arkode_mem, realtype bias)

   Specifies the bias to be applied to the error estimates within
   accuracy-based adaptivity strategies.

   **Arguments:**
      * *arkode_mem* -- pointer to the ERKStep memory block.
      * *bias* -- bias applied to error in accuracy-based time
        step estimation (default is 1.5).

   **Return value:**
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL* if the ERKStep memory is ``NULL``
      * *ARK_ILL_INPUT* if an argument has an illegal value

   **Notes:**
      Any value below 1.0 will imply a reset to the default value.



.. c:function:: int ERKStepSetFixedStepBounds(void* arkode_mem, realtype lb, realtype ub)

   Specifies the step growth interval in which the step size will remain unchanged.

   **Arguments:**
      * *arkode_mem* -- pointer to the ERKStep memory block.
      * *lb* -- lower bound on window to leave step size fixed (default is 1.0).
      * *ub* -- upper bound on window to leave step size fixed (default is 1.5).

   **Return value:**
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL* if the ERKStep memory is ``NULL``
      * *ARK_ILL_INPUT* if an argument has an illegal value

   **Notes:**
      Any interval *not* containing 1.0 will imply a reset to the default values.



.. c:function:: int ERKStepSetMaxEFailGrowth(void* arkode_mem, realtype etamxf)

   Specifies the maximum step size growth factor upon multiple successive
   accuracy-based error failures in the solver.

   **Arguments:**
      * *arkode_mem* -- pointer to the ERKStep memory block.
      * *etamxf* -- time step reduction factor on multiple error fails (default is 0.3).

   **Return value:**
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL* if the ERKStep memory is ``NULL``
      * *ARK_ILL_INPUT* if an argument has an illegal value

   **Notes:**
      Any value outside the interval :math:`(0,1]` will imply a reset to the default value.



.. c:function:: int ERKStepSetMaxFirstGrowth(void* arkode_mem, realtype etamx1)

   Specifies the maximum allowed growth factor in step size following the very
   first integration step.

   **Arguments:**
      * *arkode_mem* -- pointer to the ERKStep memory block.
      * *etamx1* -- maximum allowed growth factor after the first time
        step (default is 10000.0).

   **Return value:**
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL* if the ERKStep memory is ``NULL``
      * *ARK_ILL_INPUT* if an argument has an illegal value

   **Notes:**
      Any value :math:`\le 1.0` will imply a reset to the default value.



.. c:function:: int ERKStepSetMaxGrowth(void* arkode_mem, realtype mx_growth)

   Specifies the maximum allowed growth factor in step size between
   consecutive steps in the integration process.

   **Arguments:**
      * *arkode_mem* -- pointer to the ERKStep memory block.
      * *mx_growth* -- maximum allowed growth factor between consecutive time steps (default is 20.0).

   **Return value:**
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL* if the ERKStep memory is ``NULL``
      * *ARK_ILL_INPUT* if an argument has an illegal value

   **Notes:**
      Any value :math:`\le 1.0` will imply a reset to the default
      value.



.. c:function:: int ERKStepSetMinReduction(void* arkode_mem, realtype eta_min)

   Specifies the minimum allowed reduction factor in step size between
   step attempts, resulting from a temporal error failure in the integration
   process.

   **Arguments:**
      * *arkode_mem* -- pointer to the ERKStep memory block.
      * *eta_min* -- minimum allowed reduction factor time step after an error
        test failure (default is 0.1).

   **Return value:**
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL* if the ERKStep memory is ``NULL``
      * *ARK_ILL_INPUT* if an argument has an illegal value

   **Notes:**
      Any value :math:`\ge 1.0` or :math:`\le 0.0` will imply a reset to
      the default value.



.. c:function:: int ERKStepSetSafetyFactor(void* arkode_mem, realtype safety)

   Specifies the safety factor to be applied to the accuracy-based
   estimated step.

   **Arguments:**
      * *arkode_mem* -- pointer to the ERKStep memory block.
      * *safety* -- safety factor applied to accuracy-based time step (default is 0.96).

   **Return value:**
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL* if the ERKStep memory is ``NULL``
      * *ARK_ILL_INPUT* if an argument has an illegal value

   **Notes:**
      Any non-positive parameter will imply a reset to the default
      value.



.. c:function:: int ERKStepSetSmallNumEFails(void* arkode_mem, int small_nef)

   Specifies the threshold for "multiple" successive error failures
   before the *etamxf* parameter from
   :c:func:`ERKStepSetMaxEFailGrowth()` is applied.

   **Arguments:**
      * *arkode_mem* -- pointer to the ERKStep memory block.
      * *small_nef* -- bound to determine "multiple" for *etamxf* (default is 2).

   **Return value:**
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL* if the ERKStep memory is ``NULL``
      * *ARK_ILL_INPUT* if an argument has an illegal value

   **Notes:**
      Any non-positive parameter will imply a reset to the default value.



.. c:function:: int ERKStepSetStabilityFn(void* arkode_mem, ARKExpStabFn EStab, void* estab_data)

   Sets the problem-dependent function to estimate a stable
   time step size for the explicit portion of the ODE system.

   **Arguments:**
      * *arkode_mem* -- pointer to the ERKStep memory block.
      * *EStab* -- name of user-supplied stability function.
      * *estab_data* -- pointer to user data passed to *EStab* every time
        it is called.

   **Return value:**
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL* if the ERKStep memory is ``NULL``
      * *ARK_ILL_INPUT* if an argument has an illegal value

   **Notes:**
      This function should return an estimate of the absolute
      value of the maximum stable time step for the the ODE system.  It
      is not required, since accuracy-based adaptivity may be sufficient
      for retaining stability, but this can be quite useful for problems
      where the right-hand side function :math:`f(t,y)` contains stiff
      terms.




.. _ARKODE.Usage.ERKStep.ERKStepRootfindingInput:


Rootfinding optional input functions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The following functions can be called to set optional inputs to
control the rootfinding algorithm, the mathematics of which are
described in :numref:`ARKODE.Mathematics.Rootfinding`.


.. _ARKODE.Usage.ERKStep.ERKStepRootfindingInputTable:
.. table:: Rootfinding optional input functions

   +-----------------------------------------+------------------------------------------+----------+
   | Optional input                          | Function name                            | Default  |
   +-----------------------------------------+------------------------------------------+----------+
   | Direction of zero-crossings to monitor  | :c:func:`ERKStepSetRootDirection()`      | both     |
   +-----------------------------------------+------------------------------------------+----------+
   | Disable inactive root warnings          | :c:func:`ERKStepSetNoInactiveRootWarn()` | enabled  |
   +-----------------------------------------+------------------------------------------+----------+



.. c:function:: int ERKStepSetRootDirection(void* arkode_mem, int* rootdir)

   Specifies the direction of zero-crossings to be located and returned.

   **Arguments:**
      * *arkode_mem* -- pointer to the ERKStep memory block.
      * *rootdir* -- state array of length *nrtfn*, the number of root
        functions :math:`g_i`  (the value of *nrtfn* was supplied in
        the call to :c:func:`ERKStepRootInit()`).  If ``rootdir[i] ==
        0`` then crossing in either direction for :math:`g_i` should be
        reported.  A value of +1 or -1 indicates that the solver
        should report only zero-crossings where :math:`g_i` is
        increasing or decreasing, respectively.

   **Return value:**
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL* if the ERKStep memory is ``NULL``
      * *ARK_ILL_INPUT* if an argument has an illegal value

   **Notes:**
      The default behavior is to monitor for both zero-crossing directions.



.. c:function:: int ERKStepSetNoInactiveRootWarn(void* arkode_mem)

   Disables issuing a warning if some root function appears
   to be identically zero at the beginning of the integration.

   **Arguments:**
      * *arkode_mem* -- pointer to the ERKStep memory block.

   **Return value:**
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL* if the ERKStep memory is ``NULL``

   **Notes:**
      ERKStep will not report the initial conditions as a
      possible zero-crossing (assuming that one or more components
      :math:`g_i` are zero at the initial time).  However, if it appears
      that some :math:`g_i` is identically zero at the initial time
      (i.e., :math:`g_i` is zero at the initial time *and* after the
      first step), ERKStep will issue a warning which can be disabled with
      this optional input function.





.. _ARKODE.Usage.ERKStep.InterpolatedOutput:

Interpolated output function
--------------------------------

An optional function :c:func:`ERKStepGetDky()` is available to obtain
additional values of solution-related quantities.  This function
should only be called after a successful return from
:c:func:`ERKStepEvolve()`, as it provides interpolated values either of
:math:`y` or of its derivatives (up to the 5th derivative)
interpolated to any value of :math:`t` in the last internal step taken
by :c:func:`ERKStepEvolve()`.  Internally, this "dense output" or
"continuous extension" algorithm is identical to the algorithm used for
the maximum order implicit predictors, described in
:numref:`ARKODE.Mathematics.Predictors.Max`, except that
derivatives of the polynomial model may be evaluated upon request.



.. c:function:: int ERKStepGetDky(void* arkode_mem, realtype t, int k, N_Vector dky)

   Computes the *k*-th derivative of the function
   :math:`y` at the time *t*,
   i.e., :math:`y^{(k)}(t)`, for values of the
   independent variable satisfying :math:`t_n-h_n \le t \le t_n`, with
   :math:`t_n` as current internal time reached, and :math:`h_n` is
   the last internal step size successfully used by the solver.  This
   routine uses an interpolating polynomial of degree *min(degree, 5)*,
   where *degree* is the argument provided to
   :c:func:`ERKStepSetInterpolantDegree()`.  The user may request *k* in the
   range {0,..., *min(degree, kmax)*} where *kmax* depends on the choice of
   interpolation module. For Hermite interpolants *kmax = 5* and for Lagrange
   interpolants *kmax = 3*.

   **Arguments:**
      * *arkode_mem* -- pointer to the ERKStep memory block.
      * *t* -- the value of the independent variable at which the
        derivative is to be evaluated.
      * *k* -- the derivative order requested.
      * *dky* -- output vector (must be allocated by the user).

   **Return value:**
      * *ARK_SUCCESS* if successful
      * *ARK_BAD_K* if *k* is not in the range {0,..., *min(degree, kmax)*}.
      * *ARK_BAD_T* if *t* is not in the interval :math:`[t_n-h_n, t_n]`
      * *ARK_BAD_DKY* if the *dky* vector was ``NULL``
      * *ARK_MEM_NULL* if the ERKStep memory is ``NULL``

   **Notes:**
      It is only legal to call this function after a successful
      return from :c:func:`ERKStepEvolve()`.

      A user may access the values :math:`t_n` and :math:`h_n` via the
      functions :c:func:`ERKStepGetCurrentTime()` and
      :c:func:`ERKStepGetLastStep()`, respectively.




.. _ARKODE.Usage.ERKStep.OptionalOutputs:

Optional output functions
------------------------------

ERKStep provides an extensive set of functions that can be used to
obtain solver performance information.  We organize these into groups:

#. General ERKStep output routines are in
   :numref:`ARKODE.Usage.ERKStep.ERKStepMainOutputs`,

#. Output routines regarding root-finding results are in
   :numref:`ARKODE.Usage.ERKStep.ERKStepRootOutputs`,

#. General usability routines (e.g. to print the current ERKStep
   parameters, or output the current Butcher table) are in
   :numref:`ARKODE.Usage.ERKStep.ERKStepExtraOutputs`.

Following each table, we elaborate on each function.

Some of the optional outputs, especially the various counters, can be
very useful in determining the efficiency of various methods inside
ERKStep.  For example:

* The counters *nsteps* and *nf_evals* provide a rough measure of the
  overall cost of a given run, and can be compared between runs with
  different solver options to suggest which set of options is the most
  efficient.

* The ratio *nsteps/step_attempts* can measure the quality of the
  time step adaptivity algorithm, since a poor algorithm will result
  in more failed steps, and hence a lower ratio.

It is therefore recommended that users retrieve and output these
statistics following each run, and take some time to investigate
alternate solver options that will be more optimal for their
particular problem of interest.



.. _ARKODE.Usage.ERKStep.ERKStepMainOutputs:

Main solver optional output functions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. _ARKODE.Usage.ERKStep.ERKStepMainOutputsTable:
.. table:: Main solver optional output functions

   +------------------------------------------------------+-------------------------------------------+
   | Optional output                                      | Function name                             |
   +------------------------------------------------------+-------------------------------------------+
   | Size of ERKStep real and integer workspaces          | :c:func:`ERKStepGetWorkSpace()`           |
   +------------------------------------------------------+-------------------------------------------+
   | Cumulative number of internal steps                  | :c:func:`ERKStepGetNumSteps()`            |
   +------------------------------------------------------+-------------------------------------------+
   | Actual initial time step size used                   | :c:func:`ERKStepGetActualInitStep()`      |
   +------------------------------------------------------+-------------------------------------------+
   | Step size used for the last successful step          | :c:func:`ERKStepGetLastStep()`            |
   +------------------------------------------------------+-------------------------------------------+
   | Step size to be attempted on the next step           | :c:func:`ERKStepGetCurrentStep()`         |
   +------------------------------------------------------+-------------------------------------------+
   | Current internal time reached by the solver          | :c:func:`ERKStepGetCurrentTime()`         |
   +------------------------------------------------------+-------------------------------------------+
   | Suggested factor for tolerance scaling               | :c:func:`ERKStepGetTolScaleFactor()`      |
   +------------------------------------------------------+-------------------------------------------+
   | Error weight vector for state variables              | :c:func:`ERKStepGetErrWeights()`          |
   +------------------------------------------------------+-------------------------------------------+
   | Single accessor to many statistics at once           | :c:func:`ERKStepGetStepStats()`           |
   +------------------------------------------------------+-------------------------------------------+
   | Name of constant associated with a return flag       | :c:func:`ERKStepGetReturnFlagName()`      |
   +------------------------------------------------------+-------------------------------------------+
   | No. of explicit stability-limited steps              | :c:func:`ERKStepGetNumExpSteps()`         |
   +------------------------------------------------------+-------------------------------------------+
   | No. of accuracy-limited steps                        | :c:func:`ERKStepGetNumAccSteps()`         |
   +------------------------------------------------------+-------------------------------------------+
   | No. of attempted steps                               | :c:func:`ERKStepGetNumStepAttempts()`     |
   +------------------------------------------------------+-------------------------------------------+
   | No. of calls to *f* function                         | :c:func:`ERKStepGetNumRhsEvals()`         |
   +------------------------------------------------------+-------------------------------------------+
   | No. of local error test failures that have occurred  | :c:func:`ERKStepGetNumErrTestFails()`     |
   +------------------------------------------------------+-------------------------------------------+
   | Current ERK Butcher table                            | :c:func:`ERKStepGetCurrentButcherTable()` |
   +------------------------------------------------------+-------------------------------------------+
   | Estimated local truncation error vector              | :c:func:`ERKStepGetEstLocalErrors()`      |
   +------------------------------------------------------+-------------------------------------------+
   | Single accessor to many statistics at once           | :c:func:`ERKStepGetTimestepperStats()`    |
   +------------------------------------------------------+-------------------------------------------+
   | Number of constraint test failures                   | :c:func:`ERKStepGetNumConstrFails()`      |
   +------------------------------------------------------+-------------------------------------------+




.. c:function:: int ERKStepGetWorkSpace(void* arkode_mem, long int* lenrw, long int* leniw)

   Returns the ERKStep real and integer workspace sizes.

   **Arguments:**
      * *arkode_mem* -- pointer to the ERKStep memory block.
      * *lenrw* -- the number of ``realtype`` values in the ERKStep workspace.
      * *leniw* -- the number of integer values in the ERKStep workspace.

   **Return value:**
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL* if the ERKStep memory was ``NULL``



.. c:function:: int ERKStepGetNumSteps(void* arkode_mem, long int* nsteps)

   Returns the cumulative number of internal steps taken by
   the solver (so far).

   **Arguments:**
      * *arkode_mem* -- pointer to the ERKStep memory block.
      * *nsteps* -- number of steps taken in the solver.

   **Return value:**
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL* if the ERKStep memory was ``NULL``



.. c:function:: int ERKStepGetActualInitStep(void* arkode_mem, realtype* hinused)

   Returns the value of the integration step size used on the first step.

   **Arguments:**
      * *arkode_mem* -- pointer to the ERKStep memory block.
      * *hinused* -- actual value of initial step size.

   **Return value:**
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL* if the ERKStep memory was ``NULL``

   **Notes:**
      Even if the value of the initial integration step was
      specified by the user through a call to
      :c:func:`ERKStepSetInitStep()`, this value may have been changed by
      ERKStep to ensure that the step size fell within the prescribed
      bounds :math:`(h_{min} \le h_0 \le h_{max})`, or to satisfy the
      local error test condition.



.. c:function:: int ERKStepGetLastStep(void* arkode_mem, realtype* hlast)

   Returns the integration step size taken on the last successful
   internal step.

   **Arguments:**
      * *arkode_mem* -- pointer to the ERKStep memory block.
      * *hlast* -- step size taken on the last internal step.

   **Return value:**
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL* if the ERKStep memory was ``NULL``



.. c:function:: int ERKStepGetCurrentStep(void* arkode_mem, realtype* hcur)

   Returns the integration step size to be attempted on the next internal step.

   **Arguments:**
      * *arkode_mem* -- pointer to the ERKStep memory block.
      * *hcur* -- step size to be attempted on the next internal step.

   **Return value:**
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL* if the ERKStep memory was ``NULL``



.. c:function:: int ERKStepGetCurrentTime(void* arkode_mem, realtype* tcur)

   Returns the current internal time reached by the solver.

   **Arguments:**
      * *arkode_mem* -- pointer to the ERKStep memory block.
      * *tcur* -- current internal time reached.

   **Return value:**
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL* if the ERKStep memory was ``NULL``



.. c:function:: int ERKStepGetTolScaleFactor(void* arkode_mem, realtype* tolsfac)

   Returns a suggested factor by which the user's
   tolerances should be scaled when too much accuracy has been
   requested for some internal step.

   **Arguments:**
      * *arkode_mem* -- pointer to the ERKStep memory block.
      * *tolsfac* -- suggested scaling factor for user-supplied tolerances.

   **Return value:**
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL* if the ERKStep memory was ``NULL``



.. c:function:: int ERKStepGetErrWeights(void* arkode_mem, N_Vector eweight)

   Returns the current error weight vector.

   **Arguments:**
      * *arkode_mem* -- pointer to the ERKStep memory block.
      * *eweight* -- solution error weights at the current time.

   **Return value:**
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL* if the ERKStep memory was ``NULL``

   **Notes:**
      The user must allocate space for *eweight*, that will be
      filled in by this function.



.. c:function:: int ERKStepGetStepStats(void* arkode_mem, long int* nsteps, realtype* hinused, realtype* hlast, realtype* hcur, realtype* tcur)

   Returns many of the most useful optional outputs in a single call.

   **Arguments:**
      * *arkode_mem* -- pointer to the ERKStep memory block.
      * *nsteps* -- number of steps taken in the solver.
      * *hinused* -- actual value of initial step size.
      * *hlast* -- step size taken on the last internal step.
      * *hcur* -- step size to be attempted on the next internal step.
      * *tcur* -- current internal time reached.

   **Return value:**
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL* if the ERKStep memory was ``NULL``



.. c:function:: char *ERKStepGetReturnFlagName(long int flag)

   Returns the name of the ERKStep constant corresponding to *flag*.

   **Arguments:**
      * *flag* -- a return flag from an ERKStep function.

   **Return value:**
      The return value is a string containing the name of
      the corresponding constant.



.. c:function:: int ERKStepGetNumExpSteps(void* arkode_mem, long int* expsteps)

   Returns the cumulative number of stability-limited steps
   taken by the solver (so far).

   **Arguments:**
      * *arkode_mem* -- pointer to the ERKStep memory block.
      * *expsteps* -- number of stability-limited steps taken in the solver.

   **Return value:**
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL* if the ERKStep memory was ``NULL``



.. c:function:: int ERKStepGetNumAccSteps(void* arkode_mem, long int* accsteps)

   Returns the cumulative number of accuracy-limited steps
   taken by the solver (so far).

   **Arguments:**
      * *arkode_mem* -- pointer to the ERKStep memory block.
      * *accsteps* -- number of accuracy-limited steps taken in the solver.

   **Return value:**
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL* if the ERKStep memory was ``NULL``



.. c:function:: int ERKStepGetNumStepAttempts(void* arkode_mem, long int* step_attempts)

   Returns the cumulative number of steps attempted by the solver (so far).

   **Arguments:**
      * *arkode_mem* -- pointer to the ERKStep memory block.
      * *step_attempts* -- number of steps attempted by solver.

   **Return value:**
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL* if the ERKStep memory was ``NULL``



.. c:function:: int ERKStepGetNumRhsEvals(void* arkode_mem, long int* nf_evals)

   Returns the number of calls to the user's right-hand
   side function, :math:`f` (so far).

   **Arguments:**
      * *arkode_mem* -- pointer to the ERKStep memory block.
      * *nf_evals* -- number of calls to the user's :math:`f(t,y)` function.

   **Return value:**
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL* if the ERKStep memory was ``NULL``



.. c:function:: int ERKStepGetNumErrTestFails(void* arkode_mem, long int* netfails)

   Returns the number of local error test failures that
   have occurred (so far).

   **Arguments:**
      * *arkode_mem* -- pointer to the ERKStep memory block.
      * *netfails* -- number of error test failures.

   **Return value:**
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL* if the ERKStep memory was ``NULL``



.. c:function:: int ERKStepGetCurrentButcherTable(void* arkode_mem, ARKodeButcherTable *B)

   Returns the Butcher table currently in use by the solver.

   **Arguments:**
      * *arkode_mem* -- pointer to the ERKStep memory block.
      * *B* -- pointer to the Butcher table structure.

   **Return value:**
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL* if the ERKStep memory was ``NULL``

   **Notes:**
      The :c:type:`ARKodeButcherTable` data structure is defined as a
      pointer to the following C structure:

      .. code-block:: c

         typedef struct ARKodeButcherTableMem {

           int q;           /* method order of accuracy       */
           int p;           /* embedding order of accuracy    */
           int stages;      /* number of stages               */
           realtype **A;    /* Butcher table coefficients     */
           realtype *c;     /* canopy node coefficients       */
           realtype *b;     /* root node coefficients         */
           realtype *d;     /* embedding coefficients         */

         } *ARKodeButcherTable;

      For more details see :numref:`ARKodeButcherTable`.

.. c:function:: int ERKStepGetEstLocalErrors(void* arkode_mem, N_Vector ele)

   Returns the vector of estimated local truncation errors
   for the current step.

   **Arguments:**
      * *arkode_mem* -- pointer to the ERKStep memory block.
      * *ele* -- vector of estimated local truncation errors.

   **Return value:**
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL* if the ERKStep memory was ``NULL``

   **Notes:**
      The user must allocate space for *ele*, that will be
      filled in by this function.

      The values returned in *ele* are valid only after a successful call
      to :c:func:`ERKStepEvolve()` (i.e., it returned a non-negative value).

      The *ele* vector, together with the *eweight* vector from
      :c:func:`ERKStepGetErrWeights()`, can be used to determine how the
      various components of the system contributed to the estimated local
      error test.  Specifically, that error test uses the WRMS norm of a
      vector whose components are the products of the components of these
      two vectors.  Thus, for example, if there were recent error test
      failures, the components causing the failures are those with largest
      values for the products, denoted loosely as ``eweight[i]*ele[i]``.



.. c:function:: int ERKStepGetTimestepperStats(void* arkode_mem, long int* expsteps, long int* accsteps, long int* step_attempts, long int* nf_evals, long int* netfails)

   Returns many of the most useful time-stepper statistics in a single call.

   **Arguments:**
      * *arkode_mem* -- pointer to the ERKStep memory block.
      * *expsteps* -- number of stability-limited steps taken in the solver.
      * *accsteps* -- number of accuracy-limited steps taken in the solver.
      * *step_attempts* -- number of steps attempted by the solver.
      * *nf_evals* -- number of calls to the user's :math:`f(t,y)` function.
      * *netfails* -- number of error test failures.

   **Return value:**
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL* if the ERKStep memory was ``NULL``



.. c:function:: int ERKStepGetNumConstrFails(void* arkode_mem, long int* nconstrfails)

   Returns the cumulative number of constraint test failures (so far).

   **Arguments:**
      * *arkode_mem* -- pointer to the ERKStep memory block.
      * *nconstrfails* -- number of constraint test failures.

   **Return value:**
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL* if the ERKStep memory was ``NULL``



.. _ARKODE.Usage.ERKStep.ERKStepRootOutputs:

Rootfinding optional output functions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


.. _ARKODE.Usage.ERKStep.ERKStepRootOutputsTable:
.. table:: Rootfinding optional output functions

   +--------------------------------------------------+---------------------------------+
   | Optional output                                  | Function name                   |
   +--------------------------------------------------+---------------------------------+
   | Array showing roots found                        | :c:func:`ERKStepGetRootInfo()`  |
   +--------------------------------------------------+---------------------------------+
   | No. of calls to user root function               | :c:func:`ERKStepGetNumGEvals()` |
   +--------------------------------------------------+---------------------------------+



.. c:function:: int ERKStepGetRootInfo(void* arkode_mem, int* rootsfound)

   Returns an array showing which functions were found to
   have a root.

   **Arguments:**
      * *arkode_mem* -- pointer to the ERKStep memory block.
      * *rootsfound* -- array of length *nrtfn* with the indices of the
        user functions :math:`g_i` found to have a root (the value of
        *nrtfn* was supplied in the call to
        :c:func:`ERKStepRootInit()`).  For :math:`i = 0 \ldots`
        *nrtfn*-1, ``rootsfound[i]`` is nonzero if :math:`g_i` has a
        root, and 0 if not.

   **Return value:**
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL* if the ERKStep memory was ``NULL``

   **Notes:**
      The user must allocate space for *rootsfound* prior to
      calling this function.

      For the components of :math:`g_i` for which a root was found, the
      sign of ``rootsfound[i]`` indicates the direction of
      zero-crossing.  A value of +1 indicates that :math:`g_i` is
      increasing, while a value of -1 indicates a decreasing :math:`g_i`.



.. c:function:: int ERKStepGetNumGEvals(void* arkode_mem, long int* ngevals)

   Returns the cumulative number of calls made to the
   user's root function :math:`g`.

   **Arguments:**
      * *arkode_mem* -- pointer to the ERKStep memory block.
      * *ngevals* -- number of calls made to :math:`g` so far.

   **Return value:**
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL* if the ERKStep memory was ``NULL``




.. _ARKODE.Usage.ERKStep.ERKStepExtraOutputs:

General usability functions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The following optional routines may be called by a user to inquire
about existing solver parameters, to retrieve stored Butcher tables,
write the current Butcher table, or even to test a provided Butcher
table to determine its analytical order of accuracy.  While none of
these would typically be called during the course of solving an
initial value problem, these may be useful for users wishing to better
understand ERKStep and/or specific Runge--Kutta methods.


.. _ARKODE.Usage.ERKStep.ERKStepExtraOutputsTable:
.. table:: General usability functions

   +----------------------------------------+------------------------------------+
   | Optional routine                       | Function name                      |
   +----------------------------------------+------------------------------------+
   | Output all ERKStep solver parameters   | :c:func:`ERKStepWriteParameters()` |
   +----------------------------------------+------------------------------------+
   | Output the current Butcher table       | :c:func:`ERKStepWriteButcher()`    |
   +----------------------------------------+------------------------------------+




.. c:function:: int ERKStepWriteParameters(void* arkode_mem, FILE *fp)

   Outputs all ERKStep solver parameters to the provided file pointer.

   **Arguments:**
      * *arkode_mem* -- pointer to the ERKStep memory block.
      * *fp* -- pointer to use for printing the solver parameters.

   **Return value:**
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL* if the ERKStep memory was ``NULL``

   **Notes:**
      The *fp* argument can be ``stdout`` or ``stderr``, or it
      may point to a specific file created using ``fopen``.

      When run in parallel, only one process should set a non-NULL value
      for this pointer, since parameters for all processes would be
      identical.


.. c:function:: int ERKStepWriteButcher(void* arkode_mem, FILE *fp)

   Outputs the current Butcher table to the provided file pointer.

   **Arguments:**
      * *arkode_mem* -- pointer to the ERKStep memory block.
      * *fp* -- pointer to use for printing the Butcher table.

   **Return value:**
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL* if the ERKStep memory was ``NULL``

   **Notes:**
      The *fp* argument can be ``stdout`` or ``stderr``, or it
      may point to a specific file created using ``fopen``.

      When run in parallel, only one process should set a non-NULL value
      for this pointer, since tables for all processes would be
      identical.






.. _ARKODE.Usage.ERKStep.Reinitialization:

ERKStep re-initialization function
-------------------------------------

To reinitialize the ERKStep module for the solution of a new problem,
where a prior call to :c:func:`ERKStepCreate` has been made, the
user must call the function :c:func:`ERKStepReInit()`.  The new
problem must have the same size as the previous one.  This routine
retains the current settings for all ERKstep module options and
performs the same input checking and initializations that are done in
:c:func:`ERKStepCreate`, but it performs no memory allocation as is
assumes that the existing internal memory is sufficient for the new
problem.  A call to this re-initialization routine deletes the
solution history that was stored internally during the previous
integration.  Following a successful call to
:c:func:`ERKStepReInit()`, call :c:func:`ERKStepEvolve()` again for the
solution of the new problem.

The use of :c:func:`ERKStepReInit()` requires that the number of
Runge--Kutta stages, denoted by *s*, be no larger for the new problem than
for the previous problem.  This condition is automatically fulfilled
if the method order *q* is left unchanged.

One important use of the :c:func:`ERKStepReInit()` function is in the
treating of jump discontinuities in the RHS function.  Except in cases
of fairly small jumps, it is usually more efficient to stop at each
point of discontinuity and restart the integrator with a readjusted
ODE model, using a call to this routine.  To stop when the location
of the discontinuity is known, simply make that location a value of
``tout``.  To stop when the location of the discontinuity is
determined by the solution, use the rootfinding feature.  In either
case, it is critical that the RHS function *not* incorporate the
discontinuity, but rather have a smooth extension over the
discontinuity, so that the step across it (and subsequent rootfinding,
if used) can be done efficiently.  Then use a switch within the RHS
function (communicated through ``user_data``) that can be flipped
between the stopping of the integration and the restart, so that the
restarted problem uses the new values (which have jumped).  Similar
comments apply if there is to be a jump in the dependent variable
vector.


.. c:function:: int ERKStepReInit(void* arkode_mem, ARKRhsFn f, realtype t0, N_Vector y0)

   Provides required problem specifications and re-initializes the
   ERKStep time-stepper module.

   **Arguments:**
      * *arkode_mem* -- pointer to the ERKStep memory block.
      * *f* -- the name of the C function (of type :c:func:`ARKRhsFn()`)
        defining the right-hand side function in :math:`\dot{y} = f(t,y)`.
      * *t0* -- the initial value of :math:`t`.
      * *y0* -- the initial condition vector :math:`y(t_0)`.

   **Return value:**
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL*  if the ERKStep memory was ``NULL``
      * *ARK_MEM_FAIL*  if a memory allocation failed
      * *ARK_ILL_INPUT* if an argument has an illegal value.

   **Notes:**
      All previously set options are retained but may be updated by calling
      the appropriate "Set" functions.

      If an error occurred, :c:func:`ERKStepReInit()` also
      sends an error message to the error handler function.




.. _ARKODE.Usage.ERKStep.Reset:

ERKStep reset function
----------------------

To reset the ERKStep module to a particular state :math:`(t_R,y(t_R))` for the
continued solution of a problem, where a prior
call to :c:func:`ERKStepCreate` has been made, the user must call the function
:c:func:`ERKStepReset()`.  Like :c:func:`ERKStepReInit()` this routine retains
the current settings for all ERKStep module options and performs no memory
allocations but, unlike :c:func:`ERKStepReInit()`, this routine performs only a
*subset* of the input checking and initializations that are done in
:c:func:`ERKStepCreate`. In particular this routine retains all internal
counter values and the step size/error history. Following a successful call to
:c:func:`ERKStepReset()`, call :c:func:`ERKStepEvolve()` again to continue
solving the problem. By default the next call to :c:func:`ERKStepEvolve()` will
use the step size computed by ERKStep prior to calling :c:func:`ERKStepReset()`.
To set a different step size or have ERKStep estimate a new step size use
:c:func:`ERKStepSetInitStep()`.

One important use of the :c:func:`ERKStepReset()` function is in the
treating of jump discontinuities in the RHS functions.  Except in cases
of fairly small jumps, it is usually more efficient to stop at each
point of discontinuity and restart the integrator with a readjusted
ODE model, using a call to :c:func:`ERKStepReset()`.  To stop when
the location of the discontinuity is known, simply make that location
a value of ``tout``.  To stop when the location of the discontinuity
is determined by the solution, use the rootfinding feature.  In either
case, it is critical that the RHS functions *not* incorporate the
discontinuity, but rather have a smooth extension over the
discontinuity, so that the step across it (and subsequent rootfinding,
if used) can be done efficiently.  Then use a switch within the RHS
functions (communicated through ``user_data``) that can be flipped
between the stopping of the integration and the restart, so that the
restarted problem uses the new values (which have jumped).  Similar
comments apply if there is to be a jump in the dependent variable
vector.

.. c:function:: int ERKStepReset(void* arkode_mem, realtype tR, N_Vector yR)

   Resets the current ERKStep time-stepper module state to the provided
   independent variable value and dependent variable vector.

   **Arguments:**
      * *arkode_mem* -- pointer to the ERKStep memory block.
      * *tR* -- the value of the independent variable :math:`t`.
      * *yR* -- the value of the dependent variable vector :math:`y(t_R)`.

   **Return value:**
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL*  if the ERKStep memory was ``NULL``
      * *ARK_MEM_FAIL*  if a memory allocation failed
      * *ARK_ILL_INPUT* if an argument has an illegal value.

   **Notes:**
      By default the next call to :c:func:`ERKStepEvolve()` will use the step size
      computed by ERKStep prior to calling :c:func:`ERKStepReset()`. To set a
      different step size or have ERKStep estimate a new step size use
      :c:func:`ERKStepSetInitStep()`.

      All previously set options are retained but may be updated by calling the
      appropriate "Set" functions.

      If an error occurred, :c:func:`ERKStepReset()` also sends an error message to
      the error handler function.




.. _ARKODE.Usage.ERKStep.Resizing:

ERKStep system resize function
-------------------------------------

For simulations involving changes to the number of equations and
unknowns in the ODE system (e.g. when using spatially-adaptive
PDE simulations under a method-of-lines approach), the ERKStep
integrator may be "resized" between integration steps, through calls
to the :c:func:`ERKStepResize()` function. This function modifies
ERKStep's internal memory structures to use the new problem size,
without destruction of the temporal adaptivity heuristics.  It is
assumed that the dynamical time scales before and after the vector
resize will be comparable, so that all time-stepping heuristics prior
to calling :c:func:`ERKStepResize()` remain valid after the call.  If
instead the dynamics should be recomputed from scratch, the ERKStep
memory structure should be deleted with a call to
:c:func:`ERKStepFree()`, and recreated with a call to
:c:func:`ERKStepCreate`.

To aid in the vector resize operation, the user can supply a vector
resize function that will take as input a vector with the previous
size, and transform it in-place to return a corresponding vector of
the new size.  If this function (of type :c:func:`ARKVecResizeFn()`)
is not supplied (i.e., is set to ``NULL``), then all existing vectors
internal to ERKStep will be destroyed and re-cloned from the new input
vector.

In the case that the dynamical time scale should be modified slightly
from the previous time scale, an input *hscale* is allowed, that will
rescale the upcoming time step by the specified factor.  If a value
*hscale* :math:`\le 0` is specified, the default of 1.0 will be used.



.. c:function:: int ERKStepResize(void* arkode_mem, N_Vector yR, realtype hscale, realtype tR, ARKVecResizeFn resize, void* resize_data)

   Re-sizes ERKStep with a different state vector but with comparable
   dynamical time scale.

   **Arguments:**
      * *arkode_mem* -- pointer to the ERKStep memory block.
      * *yR* -- the newly-sized solution vector, holding the current
        dependent variable values :math:`y(t_R)`.
      * *hscale* -- the desired time step scaling factor (i.e. the next
        step will be of size *h\*hscale*).
      * *tR* -- the current value of the independent variable
        :math:`t_R` (this must be consistent with *yR*).
      * *resize* -- the user-supplied vector resize function (of type
        :c:func:`ARKVecResizeFn()`.
      * *resize_data* -- the user-supplied data structure to be passed
        to *resize* when modifying internal ERKStep vectors.

   **Return value:**
      * *ARK_SUCCESS* if successful
      * *ARK_MEM_NULL*  if the ERKStep memory was ``NULL``
      * *ARK_NO_MALLOC* if *arkode_mem* was not allocated.
      * *ARK_ILL_INPUT* if an argument has an illegal value.

   **Notes:**
      If an error occurred, :c:func:`ERKStepResize()` also sends an error
      message to the error handler function.

      If inequality constraint checking is enabled a call to
      :c:func:`ERKStepResize()` will disable constraint checking. A call
      to :c:func:`ERKStepSetConstraints()` is required to re-enable constraint
      checking.


Resizing the absolute tolerance array
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If using array-valued absolute tolerances, the absolute tolerance
vector will be invalid after the call to :c:func:`ERKStepResize()`, so
the new absolute tolerance vector should be re-set **following** each
call to :c:func:`ERKStepResize()` through a new call to
:c:func:`ERKStepSVtolerances()`.

If scalar-valued tolerances or a tolerance function was specified
through either :c:func:`ERKStepSStolerances()` or
:c:func:`ERKStepWFtolerances()`, then these will remain valid and no
further action is necessary.


.. note:: For an example showing usage of the similar
          :c:func:`ARKStepResize()` routine, see the supplied serial C
          example problem, ``ark_heat1D_adapt.c``.
