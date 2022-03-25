.. ----------------------------------------------------------------
   Programmer(s): David J. Gardner @ LLNL
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

.. _ARKODE.Usage.MRIStep.UserCallable:

MRIStep User-callable functions
==================================

This section describes the functions that are called by the
user to setup and then solve an IVP using the MRIStep time-stepping
module. Some of these are required; however, starting with
:numref:`ARKODE.Usage.MRIStep.OptionalInputs`, the functions listed involve
optional inputs/outputs or restarting, and those paragraphs may be
skipped for a casual use of ARKODE's MRIStep module. In any case,
refer to the preceding section, :numref:`ARKODE.Usage.MRIStep.Skeleton`,
for the correct order of these calls.

On an error, each user-callable function returns a negative value  (or
``NULL`` if the function returns a pointer) and sends an error message
to the error handler routine, which prints the message to ``stderr``
by default. However, the user can set a file as error output or can
provide her own error handler function (see
:numref:`ARKODE.Usage.MRIStep.OptionalInputs` for details).



.. _ARKODE.Usage.MRIStep.Initialization:

MRIStep initialization and deallocation functions
------------------------------------------------------


.. c:function:: void* MRIStepCreate(ARKRhsFn fse, ARKRhsFn fsi, realtype t0, N_Vector y0, MRIStepInnerStepper stepper, SUNContext sunctx)

   This function allocates and initializes memory for a problem to
   be solved using the MRIStep time-stepping module in ARKODE.

   **Arguments:**
      * *fse* -- the name of the function (of type :c:func:`ARKRhsFn()`)
        defining the explicit slow portion of the right-hand side function in
        :math:`\dot{y} = f^E(t,y) + f^I(t,y) + f^F(t,y)`.
      * *fsi* -- the name of the function (of type :c:func:`ARKRhsFn()`)
        defining the implicit slow portion of the right-hand side function in
        :math:`\dot{y} = f^E(t,y) + f^I(t,y) + f^F(t,y)`.
      * *t0* -- the initial value of :math:`t`.
      * *y0* -- the initial condition vector :math:`y(t_0)`.
      * *stepper* -- an :c:type:`MRIStepInnerStepper` for integrating the fast
        time scale.
      * *sunctx* -- the :c:type:`SUNContext` object (see :numref:`SUNDIALS.SUNContext`)

   **Return value:**
      If successful, a pointer to initialized problem memory of type ``void*``, to
      be passed to all user-facing MRIStep routines listed below.  If unsuccessful,
      a ``NULL`` pointer will be returned, and an error message will be printed to
      ``stderr``.

   **Example usage:**

      .. code-block:: C

         /* fast (inner) and slow (outer) ARKODE objects */
         void *inner_arkode_mem = NULL;
         void *outer_arkode_mem = NULL;

         /* MRIStepInnerStepper to wrap the inner (fast) ARKStep object */
         MRIStepInnerStepper stepper = NULL;

         /* create an ARKStep object, setting fast (inner) right-hand side
            functions and the initial condition */
         inner_arkode_mem = ARKStepCreate(ffe, ffi, t0, y0, sunctx);

         /* setup ARKStep */
         . . .

         /* create MRIStepInnerStepper wrapper for the ARKStep memory block */
         flag = ARKStepCreateMRIStepInnerStepper(inner_arkode_mem, &stepper);

         /* create an MRIStep object, setting the slow (outer) right-hand side
            functions and the initial condition */
         outer_arkode_mem = MRIStepCreate(fse, fsi, t0, y0, stepper, sunctx)

   **Example codes:**
      * ``examples/arkode/C_serial/ark_brusselator_mri.c``
      * ``examples/arkode/C_serial/ark_twowaycouple_mri.c``
      * ``examples/arkode/C_serial/ark_brusselator_1D_mri.c``
      * ``examples/arkode/C_serial/ark_onewaycouple_mri.c``
      * ``examples/arkode/C_serial/ark_reaction_diffusion_mri.c``
      * ``examples/arkode/C_serial/ark_kpr_mri.c``
      * ``examples/arkode/CXX_parallel/ark_diffusion_reaction_p.cpp``


.. c:function:: void MRIStepFree(void** arkode_mem)

   This function frees the problem memory *arkode_mem* created by
   :c:func:`MRIStepCreate`.

   **Arguments:**
      * *arkode_mem* -- pointer to the MRIStep memory block.

   **Return value:**  None


.. Tolerances
.. include:: user_callable/tolerances.rest


.. Linear/Nonlinear Solvers
.. include:: user_callable/algebraic_solvers.rest


.. Rootfinding initialization
.. include:: user_callable/rootfinding_init.rest



.. _ARKODE.Usage.MRIStep.Integration:

MRIStep solver function
-------------------------

This is the central step in the solution process -- the call to perform
the integration of the IVP.  The input argument *itask* specifies one of two
modes as to where MRIStep is to return a solution.  These modes are modified if
the user has set a stop time (with a call to the optional input function
:c:func:`MRIStepSetStopTime()`) or has requested rootfinding.


.. c:function:: int MRIStepEvolve(void* arkode_mem, realtype tout, N_Vector yout, realtype *tret, int itask)

   Integrates the ODE over an interval in :math:`t`.

   **Arguments:**
      * *arkode_mem* -- pointer to the MRIStep memory block.
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
        to the solution :math:`y(tout)` by interpolation (as described
        in :numref:`ARKODE.Mathematics.Interpolation`).

        The *ARK_ONE_STEP* option tells the solver to only take a
        single internal step :math:`y_{n-1} \to y_{n}` and then return
        control back to the calling program.  If this step will
        overtake *tout* then the solver will again return an
        interpolated result; otherwise it will return a copy of the
        internal solution :math:`y_{n}` in the vector *yout*.

   **Return value:**
      * *ARK_SUCCESS* if successful.
      * *ARK_ROOT_RETURN* if :c:func:`MRIStepEvolve()` succeeded, and
        found one or more roots.  If the number of root functions,
        *nrtfn*, is greater than 1, call
        :c:func:`MRIStepGetRootInfo()` to see which :math:`g_i` were
        found to have a root at (*\*tret*).
      * *ARK_TSTOP_RETURN* if :c:func:`MRIStepEvolve()` succeeded and
        returned at *tstop*.
      * *ARK_MEM_NULL* if the *arkode_mem* argument was ``NULL``.
      * *ARK_NO_MALLOC* if *arkode_mem* was not allocated.
      * *ARK_ILL_INPUT* if one of the inputs to
        :c:func:`MRIStepEvolve()` is illegal, or some other input to
        the solver was either illegal or missing.  Details will be
        provided in the error message.  Typical causes of this failure:

        (a) A component of the error weight vector became zero during
            internal time-stepping.

        (b) The linear solver initialization function (called by the
            user after calling :c:func:`ARKStepCreate`) failed to set
            the linear solver-specific *lsolve* field in
            *arkode_mem*.

        (c) A root of one of the root functions was found both at a
            point :math:`t` and also very near :math:`t`.

      * *ARK_TOO_MUCH_WORK* if the solver took *mxstep* internal steps
        but could not reach *tout*.  The default value for *mxstep* is
        *MXSTEP_DEFAULT = 500*.
      * *ARK_CONV_FAILURE* if convergence test failures occurred
        too many times (*ark_maxncf*) during one internal time step.
      * *ARK_LINIT_FAIL* if the linear solver's initialization
        function failed.
      * *ARK_LSETUP_FAIL* if the linear solver's setup routine failed in
        an unrecoverable manner.
      * *ARK_LSOLVE_FAIL* if the linear solver's solve routine failed in
        an unrecoverable manner.
      * *ARK_VECTOROP_ERR* a vector operation error occurred.
      * *ARK_INNERSTEP_FAILED* if the inner stepper returned with an
        unrecoverable error. The value returned from the inner stepper can be
        obtained with :c:func:`MRIStepGetLastInnerStepFlag()`.
      * *ARK_INVALID_TABLE* if an invalid coupling table was provided.

   **Notes:**
      The input vector *yout* can use the same memory as the
      vector *y0* of initial conditions that was passed to
      :c:func:`MRIStepCreate`.

      In *ARK_ONE_STEP* mode, *tout* is used only on the first call, and
      only to get the direction and a rough scale of the independent
      variable.

      All failure return values are negative and so testing the return argument for
      negative values will trap all :c:func:`MRIStepEvolve()` failures.

      Since interpolation may reduce the accuracy in the reported
      solution, if full method accuracy is desired the user should issue
      a call to :c:func:`MRIStepSetStopTime()` before the call to
      :c:func:`MRIStepEvolve()` to specify a fixed stop time to
      end the time step and return to the user.  Upon return from
      :c:func:`MRIStepEvolve()`, a copy of the internal solution
      :math:`y_{n}` will be returned in the vector *yout*.  Once the
      integrator returns at a *tstop* time, any future testing for
      *tstop* is disabled (and can be re-enabled only though a new call
      to :c:func:`MRIStepSetStopTime()`).

      On any error return in which one or more internal steps were taken
      by :c:func:`MRIStepEvolve()`, the returned values of *tret* and
      *yout* correspond to the farthest point reached in the integration.
      On all other error returns, *tret* and *yout* are left unchanged
      from those provided to the routine.

.. include:: user_callable/main_inputs.rest
.. include:: user_callable/methods.rest
.. user_callable/adaptivity.rest
.. include:: user_callable/implicit_solves.rest
.. include:: user_callable/rootfinding_inputs.rest
.. include:: user_callable/interpolated.rest
.. include:: user_callable/main_outputs.rest
.. include:: user_callable/rootfinding_outputs.rest
.. include:: user_callable/reinit.rest
.. include:: user_callable/reset.rest
.. include:: user_callable/resize.rest
