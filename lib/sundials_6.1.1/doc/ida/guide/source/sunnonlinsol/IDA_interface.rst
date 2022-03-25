.. ----------------------------------------------------------------
   SUNDIALS Copyright Start
   Copyright (c) 2002-2022, Lawrence Livermore National Security
   and Southern Methodist University.
   All rights reserved.

   See the top-level LICENSE and NOTICE files for details.

   SPDX-License-Identifier: BSD-3-Clause
   SUNDIALS Copyright End
   ----------------------------------------------------------------

.. _SUNNonlinSol.IDA:

IDA SUNNonlinearSolver interface
================================

As discussed in Chapter :numref:`IDA.Mathematics` each integration step requires the
(approximate) solution of the nonlinear system

.. math::
  G(y_n) = F\left(t_n, y_n, h_{n}^{-1}\sum_{i=0}^{q}\alpha_{n,i}y_{n-i}\right) = 0.

Rather than solving this system for the new state :math:`y_n` IDA reformulates
the system to solve for the correction :math:`y_{cor}` to the predicted new
state :math:`y_{pred}` and its derivative :math:`\dot{y}_{pred}` so that
:math:`y_n = y_{pred} + y_{cor}` and :math:`\dot{y}_n = \dot{y}_{pred} +
h_{n}^{-1}\, \alpha_{n,0}\, y_{cor}`. The nonlinear system rewritten in terms of
:math:`y_{cor}` is

.. math::
   G(y_{cor}) = F\left(t_n,\, y_{pred}+y_{cor},\,
   \dot{y}_{pred} + \alpha y_{cor}\right) = 0.
   :label: IDA_res_corrector

where :math:`\alpha = h_{n}^{-1}\, \alpha_{n,0}`.

The nonlinear system function provided by IDA to the nonlinear solver module
internally updates the current value of the new state and its derivative based
on the input correction vector. The updated vectors are used when calling the
DAE residual function and when setting up linear solves (e.g., for updating the
Jacobian or preconditioner).

IDA provides several advanced functions that will not be needed by most users,
but might be useful for users who choose to provide their own implementation of
the ``SUNNonlinearSolver`` API. For example, such a user might need access to
the current :math:`y` and :math:`\dot{y}` vectors to compute Jacobian data.


.. c:function:: int IDAGetCurrentCj(void *ida_mem, realtype *cj)

   The function ``IDAGetCurrentCj`` returns the scalar :math:`c_j` which is
   proportional to the inverse of the step size (:math:`\alpha` in
   :eq:`IDA_res_corrector`).

   **Arguments:**
      * ``ida_mem`` -- pointer to the IDA memory block.
      * ``cj`` -- the value of :math:`c_j`.

   **Return value:**
      * ``IDA_SUCCESS`` -- The optional output value has been successfully set.
      * ``IDA_MEM_NULL`` -- The IDA memory block is ``NULL``.


.. c:function:: int IDAGetCurrentY(void *ida_mem, N_Vector *ycur)

   The function ``IDAGetCurrentY`` returns the current y vector.

   **Arguments:**
      * ``ida_mem`` -- pointer to the IDA memory block.
      * ``y`` -- the current :math:`y` vector.

   **Return value:**
      * ``IDA_SUCCESS`` -- The optional output value has been successfully set.
      * ``IDA_MEM_NULL`` -- The IDA memory block is ``NULL``.


.. c:function:: int IDAGetCurrentYp(void *ida_mem, N_Vector *ypcur)

   The function ``IDAGetCurrentYp`` returns the current :math:`\dot{y}` vector.

   **Arguments:**
      * ``ida_mem`` -- pointer to the IDA memory block.
      * ``yp`` -- the current :math:`\dot{y}` vector.

   **Return value:**
      * ``IDA_SUCCESS`` -- The optional output value has been successfully set.
      * ``IDA_MEM_NULL`` -- The IDA memory block is ``NULL``.


.. c:function:: int IDAGetNonlinearSystemData(void *ida_mem, realtype *tcur, N_Vector *yypred, N_Vector *yppred, N_Vector *yyn, N_Vector *ypn, N_Vector *res, realtype *cj, void **user_data)

   The function ``IDAGetNonlinearSystemData`` returns all internal data required
   to construct the current nonlinear system :eq:`IDA_res_corrector`.

   **Arguments:**
      * ``ida_mem`` -- pointer to the IDA memory block.
      * ``tcur`` -- current value of the independent variable :math:`t_n`.
      * ``yypred`` -- predicted value of :math:`y_{pred}` at :math:`t_n`.
      * ``yppred`` -- predicted value of :math:`\dot{y}_{pred}` at :math:`t_n`.
      * ``yyn`` -- the vector :math:`y_n`. This vector may not be current and may
        need to be filled (see the note below).
      * ``ypn`` -- the vector :math:`\dot{y}_n`. This vector may not be current and
        may need to be filled (see the note below).
      * ``res`` -- the resiudal function evaluated at the current time and state,
        :math:`F(t_n, y_n, \dot{y}_n)`. This vector may not be current and may need
        to be filled (see the note below).
      * ``cj`` -- the scalar :math:`c_j` which is proportional to the inverse of
        the step size (:math:`\alpha` in :eq:`IDA_res_corrector`).
      * ``user_data`` -- pointer to the user-defined data structures.

   **Return value:**
      * ``IDA_SUCCESS`` -- The optional output values have been successfully set.
      * ``IDA_MEM_NULL`` -- The IDA memory block is ``NULL``.

   **Notes:**
      This routine is intended for users who wish to attach a custom
      :c:type:`SUNNonlinSolSysFn` to an existing ``SUNNonlinearSolver`` object
      (through a call to :c:func:`SUNNonlinSolSetSysFn`) or who need access to
      nonlinear system data to compute the nonlinear system fucntion as part of a
      custom ``SUNNonlinearSolver`` object.

      When supplying a custom :c:type:`SUNNonlinSolSysFn` to an existing
      ``SUNNonlinearSolver`` object, the user should call
      :c:func:`IDAGetNonlinearSystemData` inside the nonlinear system function to
      access the requisite data for evaluting the nonlinear system function of
      their choosing. Additionlly, if the ``SUNNonlinearSolver`` object (existing
      or custom) leverages the :c:type:`SUNNonlinSolLSetupFn` and/or
      :c:type:`SUNNonlinSolLSolveFn` functions supplied by IDA (through calls to
      :c:func:`SUNNonlinSolSetLSetupFn` and :c:func:`SUNNonlinSolSetLSolveFn`
      respectively) the vectors ``yyn`` and ``ypn``, and ``res`` must be filled in
      by the user's :c:type:`SUNNonlinSolSysFn` with the current state and
      corresponding evaluation of the right-hand side function respectively i.e.,

      .. math::
         \begin{aligned}
         yyn &= y_{pred} + y_{cor}, \\
         ypn &= \dot{y}_{pred} + \alpha \dot{y}_{cor}, \\
         res &= F\left(t_{n}, y_n, \dot{y}_n\right),
         \end{aligned}

      and :math:`f_n = f\left(t_{n}, y^n\right)` where :math:`y_{cor}` was the
      first argument supplied to the :c:type:`SUNNonlinSolSysFn`. If this function
      is called as part of a custom linear solver (i.e., the default
      :c:type:`SUNNonlinSolSysFn` is used) then the vectors ``yn`` and ``fn`` are
      only current when :c:func:`IDAGetNonlinearSystemData` is called after an
      evaluation of the nonlinear system function.


.. c:function:: int IDAComputeY(void *ida_mem, N_Vector ycor, N_Vector y)

   The function computes the current :math:`y(t)` vector based on the given
   correction vector from the nonlinear solver.

   **Arguments:**
      * ``ida_mem`` -- pointer to the IDA memory block.
      * ``ycor`` -- the correction.
      * ``y`` -- the output vector.

   **Return value:**
      * ``IDA_SUCCESS`` -- The optional output value has been successfully set.
      * ``IDA_MEM_NULL`` -- The IDA memory block is ``NULL``.


.. c:function:: int IDAComputeYp(void *ida_mem, N_Vector ycor, N_Vector yp)

   The function computes :math:`\dot{y}(t)`.

   **Arguments:**
      * ``ida_mem`` -- pointer to the IDA memory block.
      * ``ycor`` -- the correction.
      * ``yp`` -- the output vector array.

   **Return value:**
      * ``IDA_SUCCESS`` -- The optional output value has been successfully set.
      * ``IDA_MEM_NULL`` -- The IDA memory block is ``NULL``.
