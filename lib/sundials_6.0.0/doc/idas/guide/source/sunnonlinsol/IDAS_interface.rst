.. ----------------------------------------------------------------
   SUNDIALS Copyright Start
   Copyright (c) 2002-2021, Lawrence Livermore National Security
   and Southern Methodist University.
   All rights reserved.

   See the top-level LICENSE and NOTICE files for details.

   SPDX-License-Identifier: BSD-3-Clause
   SUNDIALS Copyright End
   ----------------------------------------------------------------

.. _SUNNonlinSol.IDAS:

IDAS SUNNonlinearSolver interface
=================================

As discussed in Chapter :numref:`IDAS.Mathematics` each integration step requires the
(approximate) solution of the nonlinear system

.. math::
  G(y_n) = F\left(t_n, y_n, h_{n}^{-1}\sum_{i=0}^{q}\alpha_{n,i}y_{n-i}\right) = 0.

Rather than solving this system for the new state :math:`y_n` IDAS reformulates
the system to solve for the correction :math:`y_{cor}` to the predicted new
state :math:`y_{pred}` and its derivative :math:`\dot{y}_{pred}` so that
:math:`y_n = y_{pred} + y_{cor}` and :math:`\dot{y}_n = \dot{y}_{pred} +
h_{n}^{-1}\, \alpha_{n,0}\, y_{cor}`. The nonlinear system rewritten in terms of
:math:`y_{cor}` is

.. math::
   G(y_{cor}) = F\left(t_n,\, y_{pred}+y_{cor},\,
   \dot{y}_{pred} + \alpha y_{cor}\right) = 0.
   :label: IDAS_res_corrector

where :math:`\alpha = h_{n}^{-1}\, \alpha_{n,0}`.

Similarly in the forward sensitivity analysis case the nonlinear system is also
reformulated in terms of the correction to the predicted sensitivities.

The nonlinear system function provided by IDAS to the nonlinear solver module
internally updates the current value of the new state and its derivative based
on the current corretion passed to the function (as well as the sensitivities).
These values are used when calling the DAE residual function and when setting up
linear solves (e.g., for updating the Jacobian or preconditioner).

IDAS provides several advanced functions that will not be needed by most users,
but might be useful for users who choose to provide their own implementation of
the ``SUNNonlinearSolver`` API. For example, such a user might need access to
the current :math:`y` and :math:`\dot{y}` vectors to compute Jacobian data.


.. c:function:: int IDAGetCurrentCj(void *ida_mem, realtype *cj)

   The function :c:func:`IDAGetCurrentCj` returns the scalar :math:`c_j` which is
   proportional to the inverse of the step size (:math:`\alpha` in
   :eq:`IDAS_DAE_Jacobian`).

   **Arguments:**
      * ``ida_mem`` -- pointer to the IDAS memory block.
      * ``cj`` -- the value of :math:`c_j`.

   **Return value:**
      * ``IDA_SUCCESS`` -- The optional output value has been successfully set.
      * ``IDA_MEM_NULL`` -- The IDAS memory block is ``NULL``.


.. c:function:: int IDAGetCurrentY(void *ida_mem, N_Vector *ycur)

   The function :c:func:`IDAGetCurrentY` returns the current y vector.

   **Arguments:**
      * ``ida_mem`` -- pointer to the IDAS memory block.
      * ``y`` -- the current :math:`y` vector.

   **Return value:**
      * ``IDA_SUCCESS`` -- The optional output value has been successfully set.
      * ``IDA_MEM_NULL`` -- The IDAS memory block is ``NULL``.


.. c:function:: int IDAGetCurrentYp(void *ida_mem, N_Vector *ypcur)

   The function :c:func:`IDAGetCurrentYp` returns the current :math:`\dot{y}`
   vector.

   **Arguments:**
      * ``ida_mem`` -- pointer to the IDAS memory block.
      * ``yp`` -- the current :math:`\dot{y}` vector.

   **Return value:**
      * ``IDA_SUCCESS`` -- The optional output value has been successfully set.
      * ``IDA_MEM_NULL`` -- The IDAS memory block is ``NULL``.


.. c:function:: int IDAGetCurrentYSens(void * ida_mem, N_Vector ** yyS)

   The function :c:func:`IDAGetCurrentYSens` returns the current sensitivity
   vector array.

   **Arguments:**
     * ``ida_mem`` -- pointer to the IDAS memory block.
     * ``yyS`` -- pointer to the vector array that is set to the array of sensitivity vectors.

   **Return value:**
     * ``IDA_SUCCESS`` -- The optional output value has been successfully set.
     * ``IDA_MEM_NULL`` -- The ``ida_mem`` pointer is ``NULL``.


.. c:function:: int IDAGetCurrentYpSens(void * ida_mem, N_Vector ** ypS)

   The function :c:func:`IDAGetCurrentYpSens` returns the derivative the
   current  sensitivity vector array.

   **Arguments:**
     * ``ida_mem`` -- pointer to the IDAS memory block.
     * ``ypS`` -- pointer to the vector array that is set to the array of sensitivity vector derivatives.

   **Return value:**
     * ``IDA_SUCCESS`` -- The optional output value has been successfully set.
     * ``IDA_MEM_NULL`` -- The ``ida_mem`` pointer is ``NULL``.


.. c:function:: int IDAGetNonlinearSystemData(void *ida_mem, realtype *tcur, N_Vector *yypred, N_Vector *yppred, N_Vector *yyn, N_Vector *ypn, N_Vector *res, realtype *cj, void **user_data)

   The function :c:func:`IDAGetNonlinearSystemData` returns all internal data
   required to construct the current nonlinear system :eq:`IDAS_res_corrector`.

   **Arguments:**
      * ``ida_mem`` -- pointer to the IDAS memory block.
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
        the step size (:math:`\alpha` in :eq:`IDAS_res_corrector`).
      * ``user_data`` -- pointer to the user-defined data structures.

   **Return value:**
      * ``IDA_SUCCESS`` -- The optional output values have been successfully set.
      * ``IDA_MEM_NULL`` -- The IDAS memory block is ``NULL``.

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
      :c:type:`SUNNonlinSolLSolveFn` functions supplied by IDAS (through calls to
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

      where :math:`y_{cor}` was the first argument supplied to the
      :c:type:`SUNNonlinSolSysFn`. If this function
      is called as part of a custom linear solver (i.e., the default
      :c:type:`SUNNonlinSolSysFn` is used) then the vectors ``yn``, ``ypn`` and ``res`` are
      only current when :c:func:`IDAGetNonlinearSystemData` is called after an
      evaluation of the nonlinear system function.


.. c:function:: int IDAGetNonlinearSystemDataSens(void * ida_mem, realtype* tcur, N_Vector** yySpred, N_Vector** ypSpred, N_Vector** yySn, N_Vector** ypSn, realtype* cj, void** user_data)

   The function :c:func:`IDAGetNonlinearSystemDataSens` returns all internal
   sensitivity data required to construct the current nonlinear system
   :eq:`IDAS_res_corrector`.

   **Arguments:**
     * ``ida_mem`` -- pointer to the IDAS memory block.
     * ``tcur`` -- current value of the independent variable :math:`t_n`.
     * ``yySpred`` -- predicted value of :math:`yS_{i,pred}` at :math:`t_n` for :math:`i = 0 \dots N_s - 1`.
     * ``ypSpred`` -- predicted value of :math:`\dot{y}S_{i,pred}` at :math:`t_n` for :math:`i = 0 \dots N_s - 1`.
     * ``yySn`` -- the vectors :math:`yS_{i,n}`. These vectors may be not current see the note below.
     * ``ypSn`` -- the vectors :math:`\dot{y}S_{i,n}`. These vectors may be not current see the note below.
     * ``cj`` -- the scalar :math:`c_j` which is proportional to the inverse of the step size :math:`\alpha` in :eq:`IDAS_DAE_Jacobian`.
     * ``user_data`` -- pointer to the user-defined data structures

   **Return value:**
     * ``IDA_SUCCESS`` -- The optional output values have been successfully set.
     * ``IDA_MEM_NULL`` -- The ``ida_mem`` pointer is ``NULL``.

   **Notes:**
      This routine is intended for users who wish to attach a custom
      :c:type:`SUNNonlinSolSysFn` to an  existing ``SUNNonlinearSolver`` object
      (through a call to  :c:func:`SUNNonlinSolSetSysFn`) or who need access to
      nonlinear system data to  compute the nonlinear system fucntion as part of
      a custom  ``SUNNonlinearSolver`` object.  When supplying a custom
      :c:type:`SUNNonlinSolSysFn` to an existing  ``SUNNonlinearSolver`` object,
      the user should call  :c:func:`IDAGetNonlinearSystemDataSens` inside the
      nonlinear system  function to access the requisite data for evaluting the
      nonlinear system  function of their choosing. Additionlly, if the the
      vectors ``yySn`` and  ``ypSn`` are provided as additional workspace and do
      not need to be filled in  by the user's :c:type:`SUNNonlinSolSysFn`.  If
      this function is called as part of a custom linear solver (i.e., the
      default :c:type:`SUNNonlinSolSysFn` is used) then the vectors ``yySn`` and
      ``ypSn`` are only current when :c:func:`IDAGetNonlinearSystemDataSens` is
      called after an evaluation of the nonlinear system function.


.. c:function:: int IDAComputeY(void *ida_mem, N_Vector ycor, N_Vector y)

   The function computes the current :math:`y(t)` vector based on the given
   correction vector from the nonlinear solver.

   **Arguments:**
      * ``ida_mem`` -- pointer to the IDAS memory block.
      * ``ycor`` -- the correction.
      * ``y`` -- the output vector.

   **Return value:**
      * ``IDA_SUCCESS`` -- The optional output value has been successfully set.
      * ``IDA_MEM_NULL`` -- The IDAS memory block is ``NULL``.


.. c:function:: int IDAComputeYp(void *ida_mem, N_Vector ycor, N_Vector yp)

   The function computes :math:`\dot{y}(t)`  based on the given correction
   vector from the nonlinear solver.

   **Arguments:**
      * ``ida_mem`` -- pointer to the IDAS memory block.
      * ``ycor`` -- the correction.
      * ``yp`` -- the output vector array.

   **Return value:**
      * ``IDA_SUCCESS`` -- The optional output value has been successfully set.
      * ``IDA_MEM_NULL`` -- The IDAS memory block is ``NULL``.


.. c:function:: int IDAComputeYSens(void * ida_mem, N_Vector * ycorS, N_Vector * yys)

   The function computes the sensitivities based on the given correction  vector
   from the nonlinear solver.

   **Arguments:**
     * ``ida_mem`` -- pointer to the IDAS memory block.
     * ``ycorS`` -- the correction.
     * ``yyS`` -- the output vector array.

   **Return value:**
     * ``IDA_SUCCESS`` -- The optional output value has been successfully set.
     * ``IDA_MEM_NULL`` -- The ``ida_mem`` pointer is ``NULL``.


.. c:function:: int IDAComputeYpSens(void * ida_mem, N_Vector * ycorS, N_Vector * ypS)

   The function computes the sensitivity derivatives based on the  given
   correction vector from the nonlinear solver.

   **Arguments:**
     * ``ida_mem`` -- pointer to the IDAS memory block.
     * ``ycorS`` -- the correction.
     * ``ypS`` -- the output vector array.

   **Return value:**
     * ``IDA_SUCCESS`` -- The optional output value has been successfully set.
     * ``IDA_MEM_NULL`` -- The ``ida_mem`` pointer is ``NULL``.
