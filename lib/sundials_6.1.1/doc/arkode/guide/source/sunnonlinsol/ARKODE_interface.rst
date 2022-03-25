.. ----------------------------------------------------------------
   Programmer(s): David J. Gardner @ LLNL
   ----------------------------------------------------------------
   SUNDIALS Copyright Start
   Copyright (c) 2002-2022, Lawrence Livermore National Security
   and Southern Methodist University.
   All rights reserved.

   See the top-level LICENSE and NOTICE files for details.

   SPDX-License-Identifier: BSD-3-Clause
   SUNDIALS Copyright End
   ----------------------------------------------------------------

.. _SUNNonlinSol.ARKODE:

ARKODE SUNNonlinearSolver interface
====================================

As discussed in :numref:`ARKODE.Mathematics` integration steps often require the
(approximate) solution of nonlinear systems. These systems can be formulated as
the rootfinding problem

.. math::
   \begin{array}{ll}
   G(z_i) \equiv z_i - \gamma f^I\left(t^I_{n,i}, z_i\right) - a_i = 0 &\qquad\text{[$M=I$]},\\
   G(z_i) \equiv M z_i - \gamma f^I\left(t^I_{n,i}, z_i\right) - a_i = 0 &\qquad\text{[$M$ static]},\\
   G(z_i) \equiv M(t^I_{n,i}) (z_i - a_i) - \gamma f^I\left(t^I_{n,i}, z_i\right) = 0 &\qquad\text{[$M$ time-dependent]},
   \end{array}

where :math:`z_i` is the i-th stage at time :math:`t_i` and :math:`a_i` is known
data that depends on the integration method.

Alternately, the nonlinear system above may be formulated as the fixed-point
problem

.. math::
   z_i = z_i - M(t^I_{n,i})^{-1} G(z_i),

where :math:`G(z_i)` is the variant of the rootfinding problem listed above, and
:math:`M(t^I_{n,i})` may equal either :math:`M` or :math:`I`, as applicable.

Rather than solving the above nonlinear systems for the stage value :math:`z_i`
directly, ARKODE modules solve for the correction :math:`z_{cor}` to the
predicted stage value :math:`z_{pred}` so that :math:`z_i = z_{pred} + z_{cor}`.
Thus these nonlinear systems rewritten in terms of :math:`z_{cor}` are

.. math::
   \begin{array}{ll}
   G(z_{cor}) \equiv z_{cor} - \gamma f^I\left(t^I_{n,i}, z_{i}\right) - \tilde{a}_i = 0 &\qquad\text{[$M=I$]},\\
   G(z_{cor}) \equiv M z_{cor} - \gamma f^I\left(t^I_{n,i}, z_{i}\right) - \tilde{a}_i = 0 &\qquad\text{[$M$ static]},\\
   G(z_{cor}) \equiv M(t^I_{n,i}) (z_{cor} - \tilde{a}_i) - \gamma f^I\left(t^I_{n,i}, z_{i}\right) = 0 &\qquad\text{[$M$ time-dependent]},
   \end{array}
   :label: ARKODE_Residual_corrector

for the rootfinding problem and

.. math::
   z_{cor} = z_{cor} - M(t^I_{n,i})^{-1} G(z_{i}),
   :label: ARKODE_FixedPt_corrector

for the fixed-point problem.

The nonlinear system functions provided by ARKODE modules to the nonlinear
solver module internally update the current value of the stage based on the
input correction vector i.e., :math:`z_i = z_{pred} + z_{cor}`. The updated
vector :math:`z_i` is used when calling the ODE right-hand side function and
when setting up linear solves (e.g., updating the Jacobian or preconditioner).

ARKODE modules also provide several advanced functions that will not be needed
by most users, but might be useful for users who choose to provide their own
SUNNonlinSol implementation for use by ARKODE. These routines provide
access to the internal integrator data required to evaluate
:eq:`ARKODE_Residual_corrector` or :eq:`ARKODE_FixedPt_corrector`.


ARKStep advanced output functions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Two notable functions were already listed in :numref:`ARKODE.Usage.ARKStep.ARKStepMainOutputs`:

* :c:func:`ARKStepGetCurrentState` -- returns the current state vector.
  When called within the computation of a step (i.e., during a nonlinear solve)
  this is the current stage state vector :math:`z_i = z_{pred} + z_{cor}`.
  Otherwise this is the current internal solution state vector :math:`y(t)`. In
  either case the corresponding stage or solution time can be obtained from
  :c:func:`ARKStepGetCurrentTime`.

* :c:func:`ARKStepGetCurrentGamma` -- returns the current value of the scalar :math:`\gamma`.


Additional advanced output functions that are provided to aid in the construction
of user-supplied SUNNonlinSol modules are as follows.

.. c:function:: int ARKStepGetCurrentMassMatrix(void* arkode_mem, SUNMatrix* M)

   Returns the current mass matrix. For a time dependent mass matrix the
   corresponding time can be obtained from :c:func:`ARKStepGetCurrentTime`.

   **Arguments:**
      * *arkode_mem* -- pointer to the ARKStep memory block.
      * *M* -- SUNMatrix pointer that will get set to the current mass matrix
        :math:`M(t)`. If a matrix-free method is used the output is ``NULL``.

   **Return value:**
      * ``ARK_SUCCESS`` if successful.
      * ``ARK_MEM_NULL`` if the ARKStep memory was ``NULL``.


.. c:function:: int ARKStepGetNonlinearSystemData(void* arkode_mem, realtype *tcur, N_Vector *zpred, N_Vector *z, N_Vector *Fi, realtype *gamma, N_Vector *sdata, void **user_data)

   Returns all internal data required to construct the current nonlinear
   implicit system :eq:`ARKODE_Residual_corrector` or :eq:`ARKODE_FixedPt_corrector`:

   **Arguments:**
      * *arkode_mem* -- pointer to the ARKStep memory block.
      * *tcur* -- value of the independent variable corresponding to implicit
        stage, :math:`t^I_{n,i}`.
      * *zpred* -- the predicted stage vector :math:`z_{pred}` at
        :math:`t^I_{n,i}`. This vector must not be changed.
      * *z* -- the stage vector :math:`z_{i}` above. This vector may be not
        current and may need to be filled (see the note below).
      * *Fi* -- the implicit function evaluated at the current time and state,
        :math:`f^I(t^I_{n,i}, z_{i})`. This vector may be not current and may
        need to be filled (see the note below).
      * *gamma* -- current :math:`\gamma` for implicit stage calculation.
      * *sdata* -- accumulated data from previous solution and stages,
        :math:`\tilde{a}_i`. This vector must not be changed.
      * *user_data* -- pointer to the user-defined data structure (as specified
        through :c:func:`ARKStepSetUserData`, or ``NULL`` otherwise)

   **Return value:**
      * ``ARK_SUCCESS`` if successful.
      * ``ARK_MEM_NULL`` if the ARKStep memory was ``NULL``.

   .. note::

      This routine is intended for users who whish to attach a custom
      :c:type:`SUNNonlinSolSysFn` to an existing ``SUNNonlinearSolver`` object
      (through a call to :c:func:`SUNNonlinSolSetSysFn`) or who need access to
      nonlinear system data to compute the nonlinear system function as part of
      a custom ``SUNNonlinearSolver`` object.

      When supplying a custom :c:type:`SUNNonlinSolSysFn` to an existing
      ``SUNNonlinearSolver`` object, the user should call
      :c:func:`ARKStepGetNonlinearSystemData()` **inside** the nonlinear system
      function to access the requisite data for evaluting the nonlinear systen
      function of their choosing. Additionlly, if the ``SUNNonlinearSolver`` object
      (existing or custom) leverages the :c:type:`SUNNonlinSolLSetupFn` and/or
      :c:type:`SUNNonlinSolLSolveFn` functions supplied by ARKStep (through
      calls to :c:func:`SUNNonlinSolSetLSetupFn()` and
      :c:func:`SUNNonlinSolSetLSolveFn()` respectively) the vectors *z* and *Fi*
      **must be filled** in by the user's :c:type:`SUNNonlinSolSysFn` with the
      current state and corresponding evaluation of the right-hand side function
      respectively i.e.,

      .. math::
         z  &= z_{pred} + z_{cor}, \\
         Fi &= f^I\left(t^I_{n,i}, z_i\right),

      where :math:`z_{cor}` was the first argument supplied to the
      :c:type:`SUNNonlinSolSysFn`.

      If this function is called as part of a custom linear solver (i.e., the
      default :c:type:`SUNNonlinSolSysFn` is used) then the vectors *z* and
      *Fi* are only current when :c:func:`ARKStepGetNonlinearSystemData()` is
      called after an evaluation of the nonlinear system function.


.. c:function:: int ARKStepComputeState(void* arkode_mem, N_Vector zcor, N_Vector z)

   Computes the current stage state vector using the stored prediction and the
   supplied correction from the nonlinear solver i.e.,
   :math:`z_i(t) = z_{pred} + z_{cor}`.

   **Arguments:**
      * *arkode_mem* -- pointer to the ARKStep memory block.
      * *zcor* -- the correction from the nonlinear solver.
      * *z* -- on output, the current stage state vector :math:`z_i`.

   **Return value:**
      * ``ARK_SUCCESS`` if successful.
      * ``ARK_MEM_NULL`` if the ARKStep memory was ``NULL``.



MRIStep advanced output functions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Two notable functions were already listed in :numref:`ARKODE.Usage.MRIStep.MRIStepMainOutputs`:

* :c:func:`MRIStepGetCurrentState` -- returns the current state vector. When called
  within the computation of a step (i.e., during a nonlinear solve) this is the
  current stage state vector :math:`z_i = z_{pred} + z_{cor}`. Otherwise this is
  the current internal solution state vector :math:`y(t)`. In either case the
  corresponding stage or solution time can be obtained from
  :c:func:`MRIStepGetCurrentTime()`.

* :c:func:`MRIStepGetCurrentGamma` -- returns the current value of the scalar :math:`\gamma`.


Additional advanced output functions that are provided to aid in the construction
of user-supplied SUNNonlinSol modules are as follows.


.. c:function:: int MRIStepGetNonlinearSystemData(void* arkode_mem, realtype *tcur, N_Vector *zpred, N_Vector *z, N_Vector *Fi, realtype *gamma, N_Vector *sdata, void **user_data)

   Returns all internal data required to construct the current nonlinear
   implicit system :eq:`ARKODE_Residual_corrector` or :eq:`ARKODE_FixedPt_corrector`:

   **Arguments:**
      * *arkode_mem* -- pointer to the MRIStep memory block.
      * *tcur* -- value of independent variable corresponding to slow stage
        (:math:`t^S_{n,i}` above).
      * *zpred* -- predicted nonlinear solution (:math:`z_{pred}` above). This
        vector must not be changed.
      * *z* -- stage vector (:math:`z_{i}` above). This vector may be not
        current and may need to be filled (see the note below).
      * *Fi* -- memory available for evaluating the slow implicit RHS
        (:math:`f^I(t^S_{n,i}, z_{i})` above). This vector may be
        not current and may need to be filled (see the note below).
      * *gamma* -- current :math:`\gamma` for slow stage calculation.
      * *sdata* -- accumulated data from previous solution and stages
        (:math:`\tilde{a}_i` above). This vector must not be changed.
      * *user_data* -- pointer to the user-defined data structure (as specified
        through :c:func:`MRIStepSetUserData()`, or ``NULL`` otherwise).

   **Return value:**
      * ``ARK_SUCCESS`` if successful.
      * ``ARK_MEM_NULL`` if the MRIStep memory was ``NULL``.

   .. note::

      This routine is intended for users who whish to attach a custom
      :c:type:`SUNNonlinSolSysFn` to an existing ``SUNNonlinearSolver`` object
      (through a call to :c:func:`SUNNonlinSolSetSysFn()`) or who need access to
      nonlinear system data to compute the nonlinear system function as part of
      a custom ``SUNNonlinearSolver`` object.

      When supplying a custom :c:type:`SUNNonlinSolSysFn` to an existing
      ``SUNNonlinearSolver`` object, the user should call
      :c:func:`MRIStepGetNonlinearSystemData()` **inside** the nonlinear system
      function to access the requisite data for evaluting the nonlinear systen
      function of their choosing. Additionlly, if the ``SUNNonlinearSolver`` object
      (existing or custom) leverages the :c:type:`SUNNonlinSolLSetupFn` and/or
      :c:type:`SUNNonlinSolLSolveFn` functions supplied by MRIStep (through
      calls to :c:func:`SUNNonlinSolSetLSetupFn()` and
      :c:func:`SUNNonlinSolSetLSolveFn()` respectively) the vectors *z* and *F*
      **must be filled** in by the user's :c:type:`SUNNonlinSolSysFn` with the
      current state and corresponding evaluation of the right-hand side function
      respectively i.e.,

      .. math::
         z &= z_{pred} + z_{cor}, \\
         Fi &= f^I\left(t^S_{n,i}, z_i\right),

      where :math:`z_{cor}` was the first argument supplied to the
      :c:type:`SUNNonlinSolSysFn`.

      If this function is called as part of a custom linear solver (i.e., the
      default :c:type:`SUNNonlinSolSysFn` is used) then the vectors *z* and
      *Fi* are only current when :c:func:`MRIStepGetNonlinearSystemData()` is
      called after an evaluation of the nonlinear system function.


.. c:function:: int MRIStepComputeState(void* arkode_mem, N_Vector zcor, N_Vector z)

   Computes the current stage state vector using the stored prediction and the
   supplied correction from the nonlinear solver i.e.,
   :math:`z_i = z_{pred} + z_{cor}`.

   **Arguments:**
      * *arkode_mem* -- pointer to the MRIStep memory block.
      * *zcor* -- the correction from the nonlinear solver.
      * *z* -- on output, the current stage state vector :math:`z_i`.

   **Return value:**
      * ``ARK_SUCCESS`` if successful.
      * ``ARK_MEM_NULL`` if the MRIStep memory was ``NULL``.
