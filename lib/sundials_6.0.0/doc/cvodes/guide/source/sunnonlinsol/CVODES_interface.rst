.. ----------------------------------------------------------------
   SUNDIALS Copyright Start
   Copyright (c) 2002-2021, Lawrence Livermore National Security
   and Southern Methodist University.
   All rights reserved.

   See the top-level LICENSE and NOTICE files for details.

   SPDX-License-Identifier: BSD-3-Clause
   SUNDIALS Copyright End
   ----------------------------------------------------------------

.. _SUNNonlinSol.CVODES:

CVODES SUNNonlinearSolver interface
===================================

As discussed in :numref:`CVODES.Mathematics` each integration step requires the
(approximate) solution of a nonlinear system. This system can be formulated as
the rootfinding problem

.. math::
   F(y^n) \equiv y^n - h_n \beta_{n,0} f(t_n,y^n) - a_n = 0 \, ,

or as the fixed-point problem

.. math::
   G(y^n) \equiv h_n \beta_{n,0} f(t_n,y^n) + a_n = y^n \, ,

where :math:`\displaystyle a_n\equiv\sum_{i>0}(\alpha_{n,i}y^{n-i}+h_n\beta_{n,i} {\dot{y}}^{n-i})`.

Rather than solving the above nonlinear systems for the new state :math:`y^n`
CVODES reformulates the above problems to solve for the correction :math:`y_{cor}`
to the predicted new state :math:`y_{pred}` so that :math:`y^n = y_{pred} + y_{cor}`.
The nonlinear systems rewritten in terms of :math:`y_{cor}` are

.. math::
   F(y_{cor}) \equiv y_{cor} - \gamma f(t_n, y^n) - \tilde{a}_n = 0 \,
   :label: CVODES_res_corrector

for the rootfinding problem and

.. math::
   G(y_{cor}) \equiv \gamma f(t_n, y^n) + \tilde{a}_n = y_{cor} \,
   :label: CVODES_fp_corrector

for the fixed-point problem. Similarly in the forward sensitivity analysis case
the combined state and sensitivity nonlinear systems are also reformulated in
terms of the correction to the predicted state and sensitivities.

The nonlinear system functions provided by CVODES to the nonlinear solver module
internally update the current value of the new state (and the sensitvities)
based on the input correction vector(s) i.e., :math:`y^n = y_{pred} + y_{cor}`
and :math:`s_i^n = s_{i,pred} + s_{i,cor}`. The updated vector(s) are used when
calling the right-hand side function and when setting up linear solves (e.g.,
updating the Jacobian or preconditioner).

CVODES provides several advanced functions that will not be needed by most
users, but might be useful for users who choose to provide their own
implementation of the ``SUNNonlinearSolver`` API. For example, such a user might
need access to the current value of :math:`\gamma` to compute Jacobian data.

.. c:function:: int CVodeGetCurrentGamma(void* cvode_mem, realtype* gamma)

   The function :c:func:`CVodeGetCurrentGamma` returns the current value of the
   scalar :math:`\gamma`.

   **Arguments:**
     * ``cvode_mem`` -- pointer to the CVODES memory block.
     * ``gamma`` -- the current value of the scalar :math:`\gamma` appearing in the Newton equation :math:`M = I - \gamma J`.

   **Return value:**
     * ``CV_SUCCESS`` -- The optional output value has been successfully set.
     * ``CV_MEM_NULL`` -- The CVODES memory block was ``NULL``


.. c:function:: int CVodeGetCurrentState(void* cvode_mem, N_Vector* y)

   The function :c:func:`CVodeGetCurrentState` returns the current state vector. When
   called within the computation of a step (i.e., during a nonlinear solve) this
   is :math:`y^n = y_{pred} + y_{cor}`. Otherwise this is the current internal
   solution  vector :math:`y(t)`. In either case the corresponding solution time
   can be obtained  from :c:func:`CVodeGetCurrentTime`.

   **Arguments:**
     * ``cvode_mem`` -- pointer to the CVODES memory block.
     * ``y`` -- pointer that is set to the current state vector.

   **Return value:**
     * ``CV_SUCCESS`` -- The optional output value has been successfully set.
     * ``CV_MEM_NULL`` -- The CVODES memory block was ``NULL``.


.. c:function:: int CVodeGetNonlinearSystemData(void *cvode_mem, realtype *tcur, N_Vector *ypred, N_Vector *yn, N_Vector *fn, realtype *gamma, realtype *rl1, N_Vector *zn1, void **user_data)

   The function :c:func:`CVodeGetNonlinearSystemData` returns all internal data
   required to construct the current nonlinear system :eq:`CVODES_res_corrector` or
   :eq:`CVODES_fp_corrector`.

   **Arguments:**
     * ``cvode_mem`` -- pointer to the CVODES memory block.
     * ``tn`` -- current value of the independent variable :math:`t_n`.
     * ``ypred`` -- predicted state vector :math:`y_{pred}` at :math:`t_n`.
     * ``yn`` -- state vector :math:`y^n`.  This vector may be
       not current and may need to be filled (see the note below).
     * ``fn`` -- the right-hand side function evaluated at the current time
       and state, :math:`f(t_n, y^n)`. This vector may be not current
       and may need to be filled (see the note below).
     * ``gamma`` -- current value of :math:`\gamma`.
     * ``rl1`` -- a scaling factor used to compute :math:`\tilde{a}_n = \texttt{rl1 * zn1}`.
     * ``zn1`` -- a vector used to compute :math:`\tilde{a}_n = \texttt{rl1 * zn1}`.
     * ``user_data`` -- pointer to the user-defined data structures.

   **Return value:**
     * ``CV_SUCCESS`` -- The optional output values have been successfully set.
     * ``CV_MEM_NULL`` -- The CVODES memory block was ``NULL``.

   **Notes:**
      This routine is intended for users who wish to attach a custom
      :c:type:`SUNNonlinSolSysFn` to an existing ``SUNNonlinearSolver`` object
      (through a call to  :c:func:`SUNNonlinSolSetSysFn`) or who need access to
      nonlinear system data to compute the nonlinear system function as part of a
      custom  ``SUNNonlinearSolver`` object.

      When supplying a custom
      :c:type:`SUNNonlinSolSysFn` to an existing  ``SUNNonlinearSolver`` object,
      the user should call :c:func:`CVodeGetNonlinearSystemData` inside the
      nonlinear system  function to access the requisite data for evaluting
      the nonlinear system function of their choosing. Additionlly, if the
      ``SUNNonlinearSolver`` object  (existing or custom) leverages the
      :c:type:`SUNNonlinSolLSetupFn` and/or :c:type:`SUNNonlinSolLSolveFn`
      functions supplied by CVODES (through calls to :c:func:`SUNNonlinSolSetLSetupFn`
      and :c:func:`SUNNonlinSolSetLSolveFn`, respectively) the vectors ``yn``
      and ``fn`` must be filled in by the user's  :c:type:`SUNNonlinSolSysFn`
      with the current state and corresponding evaluation of the right-hand side
      function respectively i.e.,

      .. math::
         \texttt{yn} &= y_{pred} + y_{cor},\\
         \texttt{fn} &= f\left(t_{n}, y^n\right)

      where :math:`y_{cor}` was the first argument supplied to the :c:type:`SUNNonlinSolSysFn`.

      If this function is called as part of a custom linear solver (i.e., the default
      :c:type:`SUNNonlinSolSysFn` is used) then the vectors ``yn`` and ``fn``
      are only current when :c:func:`CVodeGetNonlinearSystemData` is called after
      an evaluation of the nonlinear system function.


.. c:function:: int CVodeComputeState(void* cvode_mem, N_Vector ycor, N_Vector* yn)

   The function computes the current :math:`y(t)` vector based on stored prediction
   and the given correction vector from the nonlinear solver i.e.,
   :math:`y^n = y_{pred} + y_{cor}`.

   **Arguments:**
     * ``cvode_mem`` -- pointer to the CVODES memory block.
     * ``ycor`` -- the correction.
     * ``yn`` -- the output vector.

   **Return value:**
     * ``CV_SUCCESS`` -- The optional output value has been successfully set.
     * ``CV_MEM_NULL`` -- The CVODES memory block was ``NULL``


.. c:function:: int CVodeGetCurrentStateSens(void * cvode_mem, N_Vector ** yS)

   The function :c:func:`CVodeGetCurrentStateSens` returns the current sensitivity  state vector array.

   **Arguments:**
     * ``cvode_mem`` -- pointer to the CVODES memory block.
     * ``yS`` -- pointer to the vector array that is set to the current sensitivity state vector array.

   **Return value:**
     * ``CV_SUCCESS`` -- The optional output value has been successfully set.
     * ``CV_MEM_NULL`` -- The ``cvode_mem`` pointer is ``NULL``.


.. c:function:: int CVodeGetCurrentSensSolveIndex(void * cvode_mem, int * index)

   The function :c:func:`CVodeGetCurrentSensSolveIndex` returns the index of the  current sensitivity solve when using the ``CV_STAGGERED1`` solver.

   **Arguments:**
     * ``cvode_mem`` -- pointer to the CVODES memory block.
     * ``index`` -- will be set to the index of the current sensitivity solve.

   **Return value:**
     * ``CV_SUCCESS`` -- The optional output value has been successfully set.
     * ``CV_MEM_NULL`` -- The ``cvode_mem`` pointer is ``NULL``.


.. c:function:: int CVodeGetNonlinearSystemDataSens()

   The function :c:func:`CVodeGetNonlinearSystemDataSens` returns all internal  sensitivity data required to construct the current nonlinear system :eq:`CVODES_res_corrector` or :eq:`CVODES_fp_corrector`.

   **Arguments:**
     * ``cvode_mem`` -- pointer to the CVODES memory block.
     * ``tn`` -- current value of the independent variable :math:`t_n`.
     * ``ySpred`` -- predicted state vectors :math:`yS_{i,pred}` at :math:`t_n` for :math:`i = 0 \dots N_s - 1`. This vector must not be changed.
     * ``ySn`` -- state vectors :math:`yS_i^n` for :math:`i = 0 \dots N_s - 1`. These vectors may be not current see the note below.
     * ``gamma`` -- current value of :math:`\gamma`.
     * ``rlS1`` -- a scaling factor used to compute :math:`\tilde{a}S_n =` ``rlS1 * znS1``.
     * ``znS1`` -- a vectors used to compute :math:`\tilde{a}S_{i,n} =` ``rlS1 * znS1``.
     * ``user_data`` -- pointer to the user-defined data structure.

   **Return value:**
     * ``CV_SUCCESS`` -- The optional output values have been successfully set.
     * ``CV_MEM_NULL`` -- The ``cvode_mem`` pointer is ``NULL``.

   **Notes:**
      This routine is intended for users who whish to attach a custom
      :c:type:`SUNNonlinSolSysFn` to an  existing ``SUNNonlinearSolver``
      object (through a call to  ``SUNNonlinSolSetSysFn``) or who need access to
      nonlinear system data to  compute the nonlinear system fucntion as part of
      a custom  ``SUNNonlinearSolver`` object.  When supplying a custom
      :c:type:`SUNNonlinSolSysFn` to an existing  ``SUNNonlinearSolver`` object, the
      user should call  :c:func:`CVodeGetNonlinearSystemDataSens` inside the nonlinear
      system  function used in the sensitivity nonlinear solve to access the
      requisite data  for evaluting the nonlinear system function of their
      choosing. This could be  the same function used for solving for the new
      state (the simultaneous  approach) or a different function (the staggered
      or stagggered1 approaches).  Additionlly, the vectors ``ySn`` are only
      provided as additional  worksapce and do not need to be filled in by the
      user's  :c:type:`SUNNonlinSolSysFn`.  If this function is called as part of a
      custom linear solver (i.e., the  default :c:type:`SUNNonlinSolSysFn` is used)
      then the vectors ``ySn``  are only current when
      :c:func:`CVodeGetNonlinearSystemDataSens` is called after an  evaluation of the
      nonlinear system function.


.. c:function:: int CVodeComputeStateSens(void * cvode_mem, N_Vector* yScor, N_Vector* ySn)

   The function computes the current sensitivity vector :math:`yS(t)` for all
   sensitivities based on stored prediction and the given correction vector from
   the nonlinear solver i.e., :math:`yS^n = yS_{pred} + yS_{cor}`.

   **Arguments:**
     * ``cvode_mem`` -- pointer to the CVODES memory block.
     * ``yScor`` -- the correction.
     * ``ySn`` -- the output vector.

   **Return value:**
     * ``CV_SUCCESS`` -- The optional output value has been successfully set.
     * ``CV_MEM_NULL`` -- The ``cvode_mem`` pointer is ``NULL``.


.. c:function:: int CVodeComputeStateSens1(void * cvode_mem, sunindextype idx, N_Vector yScor1, N_Vector ySn1)

   The function computes the current sensitivity vector :math:`yS_i(t)` for the
   sensitivity at the given index based on stored prediction and the given
   correction vector from the nonlinear solver i.e.,
   :math:`yS_i^n = yS_{i,pred} + yS_{i,cor}`.

   **Arguments:**
     * ``cvode_mem`` -- pointer to the CVODES memory block.
     * ``index`` -- the index of the sensitivity to update.
     * ``yScor1`` -- the correction.
     * ``ySn1`` -- the output vector.

   **Return value:**
     * ``CV_SUCCESS`` -- The optional output value has been successfully set.
     * ``CV_MEM_NULL`` -- The ``cvode_mem`` pointer is ``NULL``.
