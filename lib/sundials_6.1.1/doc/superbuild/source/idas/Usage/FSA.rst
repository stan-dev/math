.. ----------------------------------------------------------------
   SUNDIALS Copyright Start
   Copyright (c) 2002-2022, Lawrence Livermore National Security
   and Southern Methodist University.
   All rights reserved.

   See the top-level LICENSE and NOTICE files for details.

   SPDX-License-Identifier: BSD-3-Clause
   SUNDIALS Copyright End
   ----------------------------------------------------------------

.. _IDAS.Usage.FSA:

Using IDAS for Forward Sensitivity Analysis
=============================================

This chapter describes the use of IDAS to compute solution sensitivities using
forward sensitivity analysis. One of our main guiding principles was to design
the IDAS user interface for forward sensitivity analysis as an extension of
that for IVP integration. Assuming a user main program and user-defined support
routines for IVP integration have already been defined, in order to perform
forward sensitivity analysis the user only has to insert a few more calls into
the main program and (optionally) define an additional routine which computes
the residual of the sensitivity systems :eq:`IDAS_sens_eqns`. The only
departure from this philosophy is due to the :c:type:`IDAResFn` type definition.
Without changing the definition of this type, the only way to pass values of the
problem parameters to the ODE residual function is to require the user
data structure ``f_data`` to contain a pointer to the array of real parameters
:math:`p`.

IDAS uses various constants for both input and output. These are defined as
needed in this chapter, but for convenience are also listed separately in
:numref:`IDAS.Constants`.

We begin with a brief overview, in the form of a skeleton user program.
Following that are detailed descriptions of the interface to the various
user-callable routines and of the user-supplied routines that were not already
described in :numref:`IDAS.Usage.SIM` or :numref:`IDAS.Usage.Purequad`.

.. _IDAS.Usage.FSA.skeleton_sim:

A skeleton of the user's main program
-------------------------------------

The following is a skeleton of the user's main program (or calling program) as
an application of IDAS. The user program is to have these steps in the order
indicated, unless otherwise noted. For the sake of brevity, we defer many of the
details to the later sections. As in :numref:`IDAS.Usage.SIM.skeleton_sim`,
most steps are independent of the ``N_Vector``, ``SUNMatrix``,
``SUNLinearSolver``, and ``SUNNonlinearSolver`` implementations used. For the
steps that are not, refer to Chapters :numref:`NVectors`, :numref:`SUNMatrix`,
:numref:`SUNLinSol`, :numref:`SUNNonlinSol` for the specific name of the
function to be called or macro to be referenced.

First, note that no additional header files need be included for forward
sensitivity analysis beyond those for IVP solution
:numref:`IDAS.Usage.SIM.skeleton_sim`.

Differences from the user main program skeleton in
:numref:`IDAS.Usage.SIM.skeleton_sim` are bolded.

#. Initialize parallel or multi-threaded environment

#. Create the SUNDIALS context object

#. Set the vector of initial values

#. Create matrix object

#. Create linear solver object

#. Create nonlinear solver object

#. Create IDAS object

#. Initialize IDAS solver

#. Specify integration tolerances

#. Attach linear solver

#. Set linear solver optional inputs

#. Attach nonlinear solver

#. Set nonlinear solver optional inputs

#. **Initialize quadrature integration**

   If the quadrature is not sensitivity-dependent, initialize the quadrature
   integration as described in :numref:`IDAS.Usage.Purequad`. For integrating a
   problem where the quadrature depends on the forward sensitivities see
   :numref:`IDAS.Usage.FSA.quad`.

#. **Set the sensitivity initial values**

   Call :c:func:`N_VCloneVectorArray` to create ``N_Vector`` arrays ``yS0`` and
   ``ypS0`` to hold the initial values for the sensitivity vectors of :math:`y`
   and sensitivity derivative vectors of :math:`\dot{y}`, respectively.

   .. code-block:: C

      yS0  = N_VCloneVectorArray(Ns, y0);
      ypS0 = N_VCloneVectorArray(Ns, y0);

   where ``Ns`` is the number of parameters with respect to which sensitivities
   are to be computed and ``y0`` serves only to provide an ``N_Vector`` template
   for cloning.

   Then, load initial values for each sensitivity vector ``yS0[i]`` and
   sensitivity derivative vector ``ypS0[i]`` for ``i = 0,...,N_s-1``.

#. **Activate sensitivity calculations**

   Call :c:func:`IDASensInit` to activate forward sensitivity computations
   and allocate internal memory for IDAS related to sensitivity calculations.

   If a sensitivity residual function is *not* provided to
   :c:func:`IDASensInit`, then :c:func:`IDASetSensParams` *must* be called after
   :c:func:`IDASensInit` and before :c:func:`IDASolve` to provide the array of
   problem parameters with respect to which the sensitivities are computed. This
   array must also be attached to the "user data" pointer set with
   :c:func:`IDASetUserData`. Optionally, an array of scaling factors for
   difference-quotient residual computations and a mask array to select which
   parameters with respect to which the sensitivities are computed may also be
   provided to :c:func:`IDASetSensParams`.

   check :c:func:`IDASetErrFile`

#. **Set sensitivity integration tolerances (optional)**

   Call :c:func:`IDASensSStolerances` or :c:func:`IDASensSVtolerances` to set
   the sensitivity integration tolerances or :c:func:`IDASensEEtolerances` to
   have IDAS estimate tolerances for sensitivity variables based on the
   tolerances supplied for states variables.

   If sensitivity tolerances are estimated by IDAS, the results will be more
   accurate if order of magnitude is provided by setting the ``pbar`` input to
   :c:func:`IDASetSensParams`.

#. **Create sensitivity nonlinear solver**

   If using a non-default nonlinear solver (see
   :numref:`IDAS.Usage.FSA.user_callable.nonlin_solv_init`), then create the
   desired nonlinear solver object by calling the appropriate constructor
   function defined by the particular ``SUNNonlinearSolver`` implementation
   e.g.,

   .. code-block:: c

      NLSSens = SUNNonlinSol_***Sens(...);

   for the ``IDA_SIMULTANEOUS`` or ``IDA_STAGGERED`` options ``***`` is the
   name of the nonlinear solver and ``...`` are constructor specific
   arguments (see :numref:`SUNNonlinSol` for details).

#. **Attach the sensitivity nonlinear solver**

   If using a non-default nonlinear solver, then initialize the nonlinear
   solver interface by attaching the nonlinear solver object by calling
   :c:func:`IDASetNonlinearSolverSensSim` when using the ``IDA_SIMULTANEOUS``
   corrector method, :c:func:`IDASetNonlinearSolverSensStg` when using the
   ``IDA_STAGGERED`` corrector method (see
   :numref:`IDAS.Usage.FSA.user_callable.nonlin_solv_init` for details).

#. **Set sensitivity nonlinear solver optional inputs**

   Call the appropriate set functions for the selected nonlinear solver
   module to change optional inputs specific to that nonlinear solver. These
   *must* be called after :c:func:`IDASensInit` if using the default nonlinear
   solver or after attaching a new nonlinear solver to IDAS, otherwise the
   optional inputs will be overridden by IDAS defaults. See
   :numref:`SUNNonlinSol` for more information on optional inputs.

#. Specify rootfinding problem

#. **Set optional inputs**

   Call ``IDASetSens*`` routines to change from their default values any
   optional inputs that control the behavior of IDAS in computing forward
   sensitivities. See :numref:`IDAS.Usage.FSA.user_callable.optional_inputs`
   for details.

#. Correct initial values

#. Advance solution in time

#. **Extract sensitivity solution**

   After each successful return from :c:func:`IDASolve`, the solution of the
   original IVP is available in the ``y`` argument of :c:func:`IDASolve`, while
   the sensitivity solution can be extracted into ``yS`` and ``ypS`` (which
   can be the same as ``yS0`` and ``ypS0``) by calling one of the routines
   :c:func:`IDAGetSens`, :c:func:`IDAGetSens1`, :c:func:`IDAGetSensDky`, or
   :c:func:`IDAGetSensDky1`.

#. Get optional outputs

#. Deallocate memory

   Upon completion of the integration, deallocate memory for the vectors ``yS0``
   and ``yps0`` using :c:func:`N_VDestroyVectorArray`.

#. Finalize MPI, if used


.. _IDAS.Usage.FSA.user_callable:

User-callable routines for forward sensitivity analysis
-------------------------------------------------------

This section describes the IDAS functions, in addition to those presented in
:numref:`IDAS.Usage.SIM.user_callable`, that are called by the user
to setup and solve a forward sensitivity problem.

.. _IDAS.Usage.FSA.user_callable.sensi_init:

Forward sensitivity initialization and deallocation functions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Activation of forward sensitivity computation is done by calling
:c:func:`IDASensInit` or :c:func:`IDASensInit1`, depending on whether the
sensitivity residual function returns all sensitivities at once or one by
one, respectively. The form of the call to each of these routines is as follows:

.. c:function:: int IDASensInit(void * ida_mem, int Ns, int ism, IDASensResFn fS, N_Vector * yS0, N_Vector * ypS0)

   The routine :c:func:`IDASensInit` activates forward sensitivity computations and
   allocates internal memory related to sensitivity calculations.

   **Arguments:**
     * ``ida_mem`` -- pointer to the IDAS memory block returned by
       :c:func:`IDACreate`.
     * ``Ns`` -- the number of sensitivities to be computed.
     * ``ism`` -- forward sensitivity analysis!correction strategies a flag used
       to select the sensitivity solution method. Its value can be
       ``IDA_SIMULTANEOUS`` or ``IDA_STAGGERED`` :

       * In the ``IDA_SIMULTANEOUS`` approach, the state and sensitivity
         variables are corrected at the same time. If the default Newton
         nonlinear solver is used, this amounts to performing a modified Newton
         iteration on the combined nonlinear system.
       * In the ``IDA_STAGGERED`` approach, the correction step for the
         sensitivity variables takes place at the same time for all sensitivity
         equations, but only after the correction of the state variables has
         converged and the state variables have passed the local error test.

     * ``resS`` -- is the C function which computes all sensitivity ODE
       residuals at the same time. For full details see :c:type:`IDASensResFn`.
     * ``yS0`` -- a pointer to an array of ``Ns`` vectors containing the initial
       values of the sensitivities of :math:`y`.
     * ``ypS0`` -- a pointer to an array of ``Ns`` vectors containing the
       initial values of the sensitivities of :math:`\dot{y}`.

   **Return value:**
     * ``IDA_SUCCESS`` -- The call to :c:func:`IDASensInit` was successful.
     * ``IDA_MEM_NULL`` -- The IDAS memory block was not initialized through a previous call to :c:func:`IDACreate`.
     * ``IDA_MEM_FAIL`` -- A memory allocation request has failed.
     * ``IDA_ILL_INPUT`` -- An input argument to :c:func:`IDASensInit` has an illegal value.

   **Notes:**

   Passing ``fs == NULL`` indicates using the default internal difference
   quotient sensitivity residual routine and :c:func:`IDASetSensParams` *must*
   be called before :c:func:`IDASolve`.

   If an error occurred, :c:func:`IDASensInit` also sends an error message to
   the  error handler function.

In terms of the problem size :math:`N`, number of sensitivity vectors
:math:`N_s`, and maximum method order ``maxord``, the size of the real workspace
is increased as follows:

-  Base value: :math:`\texttt{lenrw} = \texttt{lenrw} + (\texttt{maxord}+5)N_s N`

-  With :c:func:`IDASensSVtolerances`: :math:`texttt{lenrw} = \texttt{lenrw} + N_s N`

the size of the integer workspace is increased as follows:

-  Base value: :math:`\texttt{leniw} = \texttt{leniw} + (\texttt{maxord}+5)N_s N_i`

-  With :c:func:`IDASensSVtolerances`: :math:`\texttt{leniw} = \texttt{leniw} + N_s N_i`

where :math:`N_i` is the number of integers in one ``N_Vector``.

The routine :c:func:`IDASensReInit`, useful during the solution of a sequence of
problems of same size, reinitializes the sensitivity-related internal memory.
The call to it must follow a call to :c:func:`IDASensInit` (and maybe a call to
:c:func:`IDAReInit`). The number ``Ns`` of sensitivities is assumed to be
unchanged since the call to the initialization function. The call to the
:c:func:`IDASensReInit` function has the form:

.. c:function:: int IDASensReInit(void * ida_mem, int ism, N_Vector * yS0, N_Vector * ypS0)

   The routine :c:func:`IDASensReInit` reinitializes forward sensitivity computations.

   **Arguments:**
     * ``ida_mem`` -- pointer to the IDAS memory block returned by :c:func:`IDACreate`.
     * ``ism`` --  forward sensitivity analysis!correction strategies a flag used to select the sensitivity solution method. Its value can be ``IDA_SIMULTANEOUS`` , ``IDA_STAGGERED`` , or ``IDA_STAGGERED1``.
     * ``yS0`` -- a pointer to an array of ``Ns`` variables of type ``N_Vector`` containing the initial values of the sensitivities.
     * ``ypS0`` -- a pointer to an array of ``Ns`` variables of type ``N_Vector`` containing the initial values of the sensitivities of :math:`\dot{y}`.

   **Return value:**
     * ``IDA_SUCCESS`` -- The call to :c:func:`IDASensReInit` was successful.
     * ``IDA_MEM_NULL`` -- The IDAS memory block was not initialized through a previous call to :c:func:`IDACreate`.
     * ``IDA_NO_SENS`` -- Memory space for sensitivity integration was not allocated through a previous call to :c:func:`IDASensInit`.
     * ``IDA_ILL_INPUT`` -- An input argument to :c:func:`IDASensReInit` has an illegal value.
     * ``IDA_MEM_FAIL`` -- A memory allocation request has failed.

   **Notes:**

   All arguments of :c:func:`IDASensReInit` are the same as those of the
   functions :c:func:`IDASensInit`.  If an error occurred,
   :c:func:`IDASensReInit` also sends a message to the error handler function.

To deallocate all forward sensitivity-related memory (allocated in a prior call
to :c:func:`IDASensInit`), the user must call

.. c:function:: void IDASensFree(void * ida_mem)

   The function :c:func:`IDASensFree` frees the memory allocated for forward
   sensitivity computations by a previous call to :c:func:`IDASensInit`.

   **Arguments:**
     * ``ida_mem`` -- pointer to the IDAS memory block returned by :c:func:`IDACreate`.

   **Return value:**
     * The function has no return value.

   **Notes:**
      In general, :c:func:`IDASensFree` need not be called by the user, as it is
      invoked automatically by :c:func:`IDAFree`.

      After a call to :c:func:`IDASensFree`, forward sensitivity computations can be
      reactivated only by calling :c:func:`IDASensInit`.


To activate and deactivate forward sensitivity calculations for successive
IDAS runs, without having to allocate and deallocate memory, the following
function is provided:

.. c:function:: int IDASensToggleOff(void * ida_mem)

   The function :c:func:`IDASensToggleOff` deactivates forward sensitivity
   calculations. It does not deallocate sensitivity-related memory.

   **Arguments:**
     * ``ida_mem`` -- pointer to the memory previously returned by :c:func:`IDACreate`.

   **Return value:**
     * ``IDA_SUCCESS`` -- :c:func:`IDASensToggleOff` was successful.
     * ``IDA_MEM_NULL`` -- ``ida_mem`` was ``NULL``.

   **Notes:**
      Since sensitivity-related memory is not deallocated, sensitivities can  be
      reactivated at a later time (using :c:func:`IDASensReInit`).


Forward sensitivity tolerance specification functions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

One of the following three functions must be called to specify the
integration tolerances for sensitivities. Note that this call must be made after
the call to :c:func:`IDASensInit`.

.. c:function:: int IDASensSStolerances(void * ida_mem, realtype reltolS, realtype* abstolS)

   The function :c:func:`IDASensSStolerances` specifies scalar relative and absolute
   tolerances.

   **Arguments:**
     * ``ida_mem`` -- pointer to the IDAS memory block returned by :c:func:`IDACreate`.
     * ``reltolS`` -- is the scalar relative error tolerance.
     * ``abstolS`` -- is a pointer to an array of length ``Ns`` containing the scalar absolute error tolerances, one for each parameter.

   **Return value:**
     * ``IDA_SUCCESS`` -- The call to :c:func:`IDASStolerances` was successful.
     * ``IDA_MEM_NULL`` -- The IDAS memory block was not initialized through a previous call to :c:func:`IDACreate`.
     * ``IDA_NO_SENS`` -- The sensitivity allocation function :c:func:`IDASensInit` has not been called.
     * ``IDA_ILL_INPUT`` -- One of the input tolerances was negative.


.. c:function:: int IDASensSVtolerances(void * ida_mem, realtype reltolS, N_Vector* abstolS)

   The function :c:func:`IDASensSVtolerances` specifies scalar relative tolerance
   and  vector absolute tolerances.

   **Arguments:**
     * ``ida_mem`` -- pointer to the IDAS memory block returned by :c:func:`IDACreate`.
     * ``reltolS`` -- is the scalar relative error tolerance.
     * ``abstolS`` -- is an array of ``Ns`` variables of type ``N_Vector``. The ``N_Vector`` from ``abstolS[is]`` specifies the vector tolerances for ``is`` -th sensitivity.

   **Return value:**
     * ``IDA_SUCCESS`` -- The call to :c:func:`IDASVtolerances` was successful.
     * ``IDA_MEM_NULL`` -- The IDAS memory block was not initialized through a previous call to :c:func:`IDACreate`.
     * ``IDA_NO_SENS`` -- The allocation function for sensitivities has not been called.
     * ``IDA_ILL_INPUT`` -- The relative error tolerance was negative or an absolute tolerance vector had a negative component.

   **Notes:**
      This choice of tolerances is important when the absolute error tolerance
      needs to  be different for each component of any vector ``yS[i]``.


.. c:function:: int IDASensEEtolerances(void * ida_mem)

   When :c:func:`IDASensEEtolerances` is called, IDAS will estimate
   tolerances for  sensitivity variables based on the tolerances supplied for
   states variables  and the scaling factors :math:`\bar p`.

   **Arguments:**
     * ``ida_mem`` -- pointer to the IDAS memory block returned by :c:func:`IDACreate`.

   **Return value:**
     * ``IDA_SUCCESS`` -- The call to :c:func:`IDASensEEtolerances` was successful.
     * ``IDA_MEM_NULL`` -- The IDAS memory block was not initialized through a previous call to :c:func:`IDACreate`.
     * ``IDA_NO_SENS`` -- The sensitivity allocation function has not been called.


.. _IDAS.Usage.FSA.user_callable.nonlin_solv_init:

Forward sensitivity nonlinear solver interface functions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

As in the pure DAE case, when computing solution sensitivities using forward
sensitivitiy analysis IDAS uses the ``SUNNonlinearSolver`` implementation of
Newton's method defined by the ``SUNNONLINSOL_NEWTON`` module (see
:numref:`SUNNonlinSol.Newton`) by default. To specify a different nonlinear
solver in IDAS, the user's program must create a ``SUNNonlinearSolver`` object
by calling the appropriate constructor routine. The user must then attach the
``SUNNonlinearSolver`` object to IDAS by calling
:c:func:`IDASetNonlinearSolverSensSim` when using the ``IDA_SIMULTANEOUS``
corrector option, or :c:func:`IDASetNonlinearSolver` and
:c:func:`IDASetNonlinearSolverSensStg` or
:c:func:`IDASetNonlinearSolverSensStg1` when using the ``IDA_STAGGERED``
as documented below.

When changing the nonlinear solver in IDAS, :c:func:`IDASetNonlinearSolver` must
be called after :c:func:`IDAInit`; similarly
:c:func:`IDASetNonlinearSolverSensSim`, :c:func:`IDASetNonlinearSolverStg`, must
be called after :c:func:`IDASensInit`. If any calls to :c:func:`IDASolve` have been
made, then IDAS will need to be reinitialized by calling :c:func:`IDAReInit` to
ensure that the nonlinear solver is initialized correctly before any subsequent
calls to :c:func:`IDASolve`.

The first argument passed to the routines
:c:func:`IDASetNonlinearSolverSensSim`, and
:c:func:`IDASetNonlinearSolverSensStg`, is the IDAS memory pointer returned by
:c:func:`IDACreate` and the second argument is the ``SUNNonlinearSolver`` object
to use for solving the nonlinear systems :eq:`IDAS_DAE_nls`. A call to this function
attaches the nonlinear solver to the main IDAS integrator.


.. c:function:: int IDASetNonlinearSolverSensSim(void * ida_mem, SUNNonlinearSolver NLS)

   The function :c:func:`IDASetNonLinearSolverSensSim` attaches a
   ``SUNNonlinearSolver``  object (``NLS``) to IDAS when using the
   ``IDA_SIMULTANEOUS`` approach to  correct the state and sensitivity variables
   at the same time.

   **Arguments:**
     * ``ida_mem`` -- pointer to the IDAS memory block.
     * ``NLS`` -- ``SUNNonlinearSolver`` object to use for solving nonlinear
       system :eq:`IDAS_DAE_nls`.

   **Return value:**
     * ``IDA_SUCCESS`` -- The nonlinear solver was successfully attached.
     * ``IDA_MEM_NULL`` -- The ``ida_mem`` pointer is ``NULL``.
     * ``IDA_ILL_INPUT`` -- The SUNNONLINSOL object is ``NULL`` , does not implement the required nonlinear solver operations, is not of the correct type, or the residual function, convergence test function, or maximum number of nonlinear iterations could not be set.


.. c:function:: int IDASetNonlinearSolverSensStg(void * ida_mem, SUNNonlinearSolver NLS)

   The function :c:func:`IDASetNonLinearSolverSensStg` attaches a
   ``SUNNonlinearSolver``  object (``NLS``) to IDAS when using the
   ``IDA_STAGGERED`` approach to  correct all the sensitivity variables after the
   correction of the state  variables.

   **Arguments:**
     * ``ida_mem`` -- pointer to the IDAS memory block.
     * ``NLS`` -- SUNNONLINSOL object to use for solving nonlinear systems.

   **Return value:**
     * ``IDA_SUCCESS`` -- The nonlinear solver was successfully attached.
     * ``IDA_MEM_NULL`` -- The ``ida_mem`` pointer is ``NULL``.
     * ``IDA_ILL_INPUT`` -- The SUNNONLINSOL object is ``NULL`` , does not implement the required nonlinear solver operations, is not of the correct type, or the residual function, convergence test function, or maximum number of nonlinear iterations could not be set.

   **Notes:**
      This function only attaches the ``SUNNonlinearSolver`` object for
      correcting the  sensitivity variables. To attach a ``SUNNonlinearSolver``
      object for the state  variable correction use
      :c:func:`IDASetNonlinearSolver`.

Forward sensitivity initial condition calculation function
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

:c:func:`IDACalcIC` also calculates corrected initial conditions for sensitivity
variables of a DAE system. When used for initial conditions calculation of the
forward sensitivities, :c:func:`IDACalcIC` must be preceded by successful calls
to :c:func:`IDASensInit` (or :c:func:`IDASensReInit`) and should precede the
call(s) to :c:func:`IDASolve`. For restrictions that apply for initial
conditions calculation of the state variables, see
:numref:`IDAS.Usage.SIM.user_callable.initialcondition`.

Calling :c:func:`IDACalcIC` is optional. It is only necessary when the initial
conditions do not satisfy the sensitivity systems. Even if forward sensitivity
analysis was enabled, the call to the initial conditions calculation function
:c:func:`IDACalcIC` is exactly the same as for state variables.

.. code-block:: C

   flag = IDACalcIC(ida_mem, icopt, tout1);

See :c:func:`IDACalcIC` for a list of possible return values.


IDAS solver function
^^^^^^^^^^^^^^^^^^^^^^

Even if forward sensitivity analysis was enabled, the call to the main solver
function :c:func:`IDASolve` is exactly the same as in :numref:`IDAS.Usage.SIM`.
However, in this case the return value ``flag`` can also be one of the
following:

- ``IDA_SRES_FAIL`` -- The sensitivity residual function failed in an
  unrecoverable manner.
- ``IDA_REP_SRES_ERR`` -- The user's residual function
  repeatedly returned a recoverable error flag, but
  the solver was unable to recover.

.. _IDAS.Usage.FSA.user_callable.sensi_get:

Forward sensitivity extraction functions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If forward sensitivity computations have been initialized by a call to
:c:func:`IDASensInit`, or reinitialized by a call to :c:func:`IDASensReInit`,
then IDAS computes both a solution and sensitivities at time ``t``. However,
:c:func:`IDASolve` will still return only the solution :math:`y` in ``yout``.
Solution sensitivities can be obtained through one of the following functions:


.. c:function:: int IDAGetSens(void * ida_mem, realtype * tret, N_Vector * yS)

   The function :c:func:`IDAGetSens` returns the sensitivity solution vectors after
   a  successful return from :c:func:`IDASolve`.

   **Arguments:**
     * ``ida_mem`` -- pointer to the memory previously allocated by :c:func:`IDAInit`.
     * ``tret`` -- the time reached by the solver output.
     * ``yS`` -- array of computed forward sensitivity vectors. This vector array must be allocated by the user.

   **Return value:**
     * ``IDA_SUCCESS`` -- :c:func:`IDAGetSens` was successful.
     * ``IDA_MEM_NULL`` -- ``ida_mem`` was ``NULL``.
     * ``IDA_NO_SENS`` -- Forward sensitivity analysis was not initialized.
     * ``IDA_BAD_DKY`` -- ``yS`` is ``NULL``.

   **Notes:**
      Note that the argument ``tret`` is an output for this function. Its value
      will be the same as that returned at the last :c:func:`IDASolve` call.


The function :c:func:`IDAGetSensDky` computes the ``k``-th derivatives of the
interpolating polynomials for the sensitivity variables at time ``t``. This
function is called by :c:func:`IDAGetSens` with ``k`` :math:`= 0`, but may also be
called directly by the user.


.. c:function:: int IDAGetSensDky(void * ida_mem, realtype t, int k, N_Vector * dkyS)

   The function :c:func:`IDAGetSensDky` returns derivatives of the sensitivity
   solution  vectors after a successful return from :c:func:`IDASolve`.

   **Arguments:**
     * ``ida_mem`` -- pointer to the memory previously allocated by :c:func:`IDAInit`.
     * ``t`` -- specifies the time at which sensitivity information is requested. The time ``t`` must fall within the interval defined by the last successful step taken by IDAS.
     * ``k`` -- order of derivatives. ``k`` must be in the range :math:`0, 1, ..., klast` where :math:`klast` is the method order of the last successful step.
     * ``dkyS`` -- array of ``Ns`` vectors containing the derivatives on output. The space for ``dkyS`` must be allocated by the user.

   **Return value:**
     * ``IDA_SUCCESS`` -- :c:func:`IDAGetSensDky` succeeded.
     * ``IDA_MEM_NULL`` -- ``ida_mem`` was ``NULL``.
     * ``IDA_NO_SENS`` -- Forward sensitivity analysis was not initialized.
     * ``IDA_BAD_DKY`` -- One of the vectors ``dkyS[i]`` is ``NULL``.
     * ``IDA_BAD_K`` -- ``k`` is not in the range :math:`0, 1, ...,` ``qlast``.
     * ``IDA_BAD_T`` -- The time ``t`` is not in the allowed range.


Forward sensitivity solution vectors can also be extracted separately for each
parameter in turn through the functions :c:func:`IDAGetSens1` and
:c:func:`IDAGetSensDky1`, defined as follows:


.. c:function:: int IDAGetSens1(void * ida_mem, realtype * tret, int is, N_Vector yS)

   The function ``IDAGetSens1`` returns the ``is``-th sensitivity solution
   vector  after a successful return from :c:func:`IDASolve`.

   **Arguments:**
     * ``ida_mem`` -- pointer to the memory previously allocated by :c:func:`IDAInit`.
     * ``tret`` -- the time reached by the solver output.
     * ``is`` -- specifies which sensitivity vector is to be returned :math:`0\le` ``is`` :math:`< N_s`.
     * ``yS`` -- the computed forward sensitivity vector. This vector array must be allocated by the user.

   **Return value:**
     * ``IDA_SUCCESS`` -- ``IDAGetSens1`` was successful.
     * ``IDA_MEM_NULL`` -- ``ida_mem`` was ``NULL``.
     * ``IDA_NO_SENS`` -- Forward sensitivity analysis was not initialized.
     * ``IDA_BAD_IS`` -- The index ``is`` is not in the allowed range.
     * ``IDA_BAD_DKY`` -- ``yS`` is ``NULL``.
     * ``IDA_BAD_T`` -- The time ``t`` is not in the allowed range.

   **Notes:**
      Note that the argument ``tret`` is an output for this function. Its value
      will be  the same as that returned at the last :c:func:`IDASolve` call.


.. c:function:: int IDAGetSensDky1(void * ida_mem, realtype t, int k, int is, N_Vector dkyS)

   The function ``IDAGetSensDky1`` returns the ``k``-th derivative of the
   ``is``-th sensitivity solution vector after a successful return from
   :c:func:`IDASolve`.

   **Arguments:**
     * ``ida_mem`` -- pointer to the memory previously allocated by :c:func:`IDAInit`.
     * ``t`` -- specifies the time at which sensitivity information is requested. The time ``t`` must fall within the interval defined by the last successful step taken by IDAS.
     * ``k`` -- order of derivative.
     * ``is`` -- specifies the sensitivity derivative vector to be returned :math:`0\le` ``is`` :math:`< N_s`.
     * ``dkyS`` -- the vector containing the derivative. The space for ``dkyS`` must be allocated by the user.

   **Return value:**
     * ``IDA_SUCCESS`` -- ``IDAGetQuadDky1`` succeeded.
     * ``IDA_MEM_NULL`` -- The pointer to ``ida_mem`` was ``NULL``.
     * ``IDA_NO_SENS`` -- Forward sensitivity analysis was not initialized.
     * ``IDA_BAD_DKY`` -- ``dkyS`` or one of the vectors ``dkyS[i]`` is ``NULL``.
     * ``IDA_BAD_IS`` -- The index ``is`` is not in the allowed range.
     * ``IDA_BAD_K`` -- ``k`` is not in the range :math:`0, 1, ...,` ``qlast``.
     * ``IDA_BAD_T`` -- The time ``t`` is not in the allowed range.


.. _IDAS.Usage.FSA.user_callable.optional_inputs:

Optional inputs for forward sensitivity analysis
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Optional input variables that control the computation of sensitivities can be
changed from their default values through calls to ``IDASetSens*`` functions.
:numref:`IDAS.Usage.FSA.user_callable.optional_inputs.Table` lists all forward
sensitivity optional input functions in IDAS which are described in detail in
the remainder of this section.

We note that, on an error return, all of the optional input functions send an
error message to the error handler function. All error return values are
negative, so the test ``flag < 0`` will catch all errors. Finally, a call to a
``IDASetSens***`` function can be made from the user's calling program at any
time and, if successful, takes effect immediately.

.. _IDAS.Usage.FSA.user_callable.optional_inputs.Table:
.. table:: Forward sensitivity optional inputs
   :align: center

  =================================== ==================================== ============
  **Optional input**                  **Routine name**                     **Default**
  =================================== ==================================== ============
  Sensitivity scaling factors         :c:func:`IDASetSensParams`           ``NULL``
  DQ approximation method             :c:func:`IDASetSensDQMethod`         centered/0.0
  Error control strategy              :c:func:`IDASetSensErrCon`           ``SUNFALSE``
  Maximum no. of nonlinear iterations :c:func:`IDASetSensMaxNonlinIters`   4
  =================================== ==================================== ============


.. c:function:: int IDASetSensParams(void * ida_mem, realtype * p, realtype * pbar, int * plist)

   The function :c:func:`IDASetSensParams` specifies problem parameter information
   for sensitivity calculations.

   **Arguments:**
     * ``ida_mem`` -- pointer to the IDAS memory block.
     * ``p`` -- a pointer to the array of real problem parameters used to
       evaluate :math:`F(t,y,\dot{y},p)`. If non- ``NULL`` , ``p`` must point to
       a field in the user's data structure ``user_data`` passed to the residual
       function.
     * ``pbar`` -- an array of ``Ns`` positive scaling factors.
       If non- ``NULL`` , ``pbar`` must have all its components :math:`> 0.0`.
     * ``plist`` -- an array of ``Ns`` non-negative indices to specify
       which components ``p[i]`` to use in estimating the sensitivity equations.
       If non- ``NULL`` , ``plist`` must have all components :math:`\ge 0`.

   **Return value:**
     * ``IDA_SUCCESS`` -- The optional value has been successfully set.
     * ``IDA_MEM_NULL`` -- The ``ida_mem`` pointer is ``NULL``.
     * ``IDA_NO_SENS`` -- Forward sensitivity analysis was not initialized.
     * ``IDA_ILL_INPUT`` -- An argument has an illegal value.

   .. note::

      The array ``p`` only needs to include the parameters with respect to which
      sensitivities are (potentially) desired.

      If the user provides a function to evaluate the sensitivity residuals,
      ``p`` need not be specified.

      When computing the sensitivity residual via a difference-quotient or
      estimating sensitivity tolerances the results will be more accurate if
      order of magnitude information is provided with ``pbar``. Typically, if
      ``p[0] != 0``, the value ``pbar[i] = abs(p[plist[i]])`` can be used. By
      default IDAS uses ``p[i] = 1.0``.

      If the user provides a function to evaluate the sensitivity residual and
      specifies tolerances for the sensitivity variables, ``pbar`` need not be
      specified.

      By default IDA computes sensitivities with respect to the first ``Ns``
      parameters in ``p`` i.e., ``plist[i] = i`` for ``i = 0,...,Ns-1``. If
      sensitivities with respect to the :math:`j`-th parameter ``p[j]`` are
      desired, set ``plist[i] = j`` for some :math:`0 \leq i < N_s` and
      :math:`0 \leq j < N_p` where :math:`N_p` is the number of element in
      ``p``.

      If the user provides a function to evaluate the sensitivity residuals,
      ``plist`` need not be specified.

   .. warning::

      This function must be preceded by a call to :c:func:`IDASensInit`.

      The array ``p`` *must* also be attached to the user data structure. For
      example, ``user_data->p = p;``.

.. c:function:: int IDASetSensDQMethod(void * ida_mem, int DQtype, realtype DQrhomax)

   The function :c:func:`IDASetSensDQMethod` specifies the difference quotient
   strategy in  the case in which the residual of the sensitivity
   equations are to  be computed by IDAS.

   **Arguments:**
     * ``ida_mem`` -- pointer to the IDAS memory block.
     * ``DQtype`` -- specifies the difference quotient type. Its value can be ``IDA_CENTERED`` or ``IDA_FORWARD``.
     * ``DQrhomax`` -- positive value of the selection parameter used in deciding switching between a simultaneous or separate approximation of the two terms in the sensitivity residual.

   **Return value:**
     * ``IDA_SUCCESS`` -- The optional value has been successfully set.
     * ``IDA_MEM_NULL`` -- The ``ida_mem`` pointer is ``NULL``.
     * ``IDA_ILL_INPUT`` -- An argument has an illegal value.

   **Notes:**

   If ``DQrhomax`` :math:`= 0.0`, then no switching is performed. The
   approximation is done simultaneously using either centered or forward finite
   differences, depending on the value of ``DQtype``.  For values of
   ``DQrhomax`` :math:`\ge 1.0`, the simultaneous approximation is used whenever
   the estimated finite difference perturbations for states and parameters are
   within a factor of ``DQrhomax``, and the separate approximation is used
   otherwise. Note that a value ``DQrhomax`` :math:`<1.0` will effectively
   disable switching.  See :numref:`IDAS.Mathematics.FSA` for more details.

   The default value are ``DQtype == IDA_CENTERED`` and
   ``DQrhomax``:math:`=0.0`.


.. c:function:: int IDASetSensErrCon(void * ida_mem, booleantype errconS)

   The function :c:func:`IDASetSensErrCon` specifies the error control  strategy for
   sensitivity variables.

   **Arguments:**
     * ``ida_mem`` -- pointer to the IDAS memory block.
     * ``errconS`` -- specifies whether sensitivity variables are to be included ``SUNTRUE`` or not ``SUNFALSE`` in the error control mechanism.

   **Return value:**
     * ``IDA_SUCCESS`` -- The optional value has been successfully set.
     * ``IDA_MEM_NULL`` -- The ``ida_mem`` pointer is ``NULL``.

   **Notes:**
      By default, ``errconS`` is set to ``SUNFALSE``.  If ``errconS = SUNTRUE``
      then both state variables and  sensitivity variables are included in the
      error tests.  If ``errconS = SUNFALSE`` then the sensitivity
      variables are excluded from the  error tests. Note that, in any event, all
      variables are considered in the convergence  tests.


.. c:function:: int IDASetSensMaxNonlinIters(void * ida_mem, int maxcorS)

   The function :c:func:`IDASetSensMaxNonlinIters` specifies the maximum  number of
   nonlinear solver iterations for sensitivity variables per step.

   **Arguments:**
     * ``ida_mem`` -- pointer to the IDAS memory block.
     * ``maxcorS`` -- maximum number of nonlinear solver iterations allowed per step :math:`> 0`.

   **Return value:**
     * ``IDA_SUCCESS`` -- The optional value has been successfully set.
     * ``IDA_MEM_NULL`` -- The ``ida_mem`` pointer is ``NULL``.
     * ``IDA_MEM_FAIL`` -- The SUNNONLINSOL module is ``NULL``.

   **Notes:**
      The default value is 3.


.. _IDAS.Usage.FSA.user_callable.optional_output:

Optional outputs for forward sensitivity analysis
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Optional output functions that return statistics and solver performance
information related to forward sensitivity computations are listed in
:numref:`IDAS.Usage.FSA.user_callable.optional_output.Table` and described in
detail in the remainder of this section.

.. _IDAS.Usage.FSA.user_callable.optional_output.Table:
.. table:: Forward sensitivity optional outputs
   :align: center

   ================================================== ================================================
   **Optional output**                                **Routine name**
   ================================================== ================================================
   No. of calls to sensitivity residual function      :c:func:`IDAGetSensNumResEvals`
   No. of calls to residual function for sensitivity  :c:func:`IDAGetNumResEvalsSens`
   No. of sensitivity local error test failures       :c:func:`IDAGetSensNumErrTestFails`
   No. of calls to lin. solv. setup routine for sens. :c:func:`IDAGetSensNumLinSolvSetups`
   Error weight vector for sensitivity variables      :c:func:`IDAGetSensErrWeights`
   Sensitivity-related statistics as a group          :c:func:`IDAGetSensStats`
   No. of sens. nonlinear solver iterations           :c:func:`IDAGetSensNumNonlinSolvIters`
   No. of sens. convergence failures                  :c:func:`IDAGetSensNumNonlinSolvConvFails`
   Sens. nonlinear solver statistics as a group       :c:func:`IDAGetSensNonlinSolveStats`
   ================================================== ================================================


.. c:function:: int IDAGetSensNumResEvals(void * ida_mem, long int * nfSevals)

   The function :c:func:`IDAGetSensNumResEvals` returns the number of calls to the
   sensitivity  residual function.

   **Arguments:**
     * ``ida_mem`` -- pointer to the IDAS memory block.
     * ``nfSevals`` -- number of calls to the sensitivity residual function.

   **Return value:**
     * ``IDA_SUCCESS`` -- The optional output value has been successfully set.
     * ``IDA_MEM_NULL`` -- The ``ida_mem`` pointer is ``NULL``.
     * ``IDA_NO_SENS`` -- Forward sensitivity analysis was not initialized.


.. c:function:: int IDAGetNumResEvalsSens(void * ida_mem, long int * nfevalsS)

   The function :c:func:`IDAGetNumResEvalsSEns` returns the number of calls to the
   user's residual function due to the internal finite difference
   approximation  of the sensitivity residuals.

   **Arguments:**
     * ``ida_mem`` -- pointer to the IDAS memory block.
     * ``nfevalsS`` -- number of calls to the user's DAE residual function for the evaluation of sensitivity residuals.

   **Return value:**
     * ``IDA_SUCCESS`` -- The optional output value has been successfully set.
     * ``IDA_MEM_NULL`` -- The ``ida_mem`` pointer is ``NULL``.
     * ``IDA_NO_SENS`` -- Forward sensitivity analysis was not initialized.

   **Notes:**
      This counter is incremented only if the internal finite difference
      approximation  routines are used for the evaluation of the sensitivity
      residuals.


.. c:function:: int IDAGetSensNumErrTestFails(void * ida_mem, long int * nSetfails)

   The function :c:func:`IDAGetSensNumErrTestFails` returns the number of local
   error test failures for the sensitivity variables that have occurred.

   **Arguments:**
     * ``ida_mem`` -- pointer to the IDAS memory block.
     * ``nSetfails`` -- number of error test failures.

   **Return value:**
     * ``IDA_SUCCESS`` -- The optional output value has been successfully set.
     * ``IDA_MEM_NULL`` -- The ``ida_mem`` pointer is ``NULL``.
     * ``IDA_NO_SENS`` -- Forward sensitivity analysis was not initialized.

   **Notes:**
      This counter is incremented only if the sensitivity variables have been
      included in the error test (see :c:func:`IDASetSensErrCon`).  Even in
      that case, this counter is not incremented if the
      ``ism = IDA_SIMULTANEOUS``  sensitivity solution method has been used.


.. c:function:: int IDAGetSensNumLinSolvSetups(void * ida_mem, long int * nlinsetupsS)

   The function :c:func:`IDAGetSensNumLinSolvSetups` returns the number of calls  to the linear solver setup function due to forward sensitivity calculations.

   **Arguments:**
     * ``ida_mem`` -- pointer to the IDAS memory block.
     * ``nlinsetupsS`` -- number of calls to the linear solver setup function.

   **Return value:**
     * ``IDA_SUCCESS`` -- The optional output value has been successfully set.
     * ``IDA_MEM_NULL`` -- The ``ida_mem`` pointer is ``NULL``.
     * ``IDA_NO_SENS`` -- Forward sensitivity analysis was not initialized.

   **Notes:**
      This counter is incremented only if a nonlinear solver requiring a linear
      solve has been used and the ``ism = IDA_STAGGERED`` sensitivity solution
      method has been specified (see
      :numref:`IDAS.Usage.FSA.user_callable.sensi_init`).


.. c:function:: int IDAGetSensStats(void* ida_mem, long int* nresSevals, \
                long int* nresevalsS, long int* nSetfails, \
                long int* nlinsetupsS)

   The function :c:func:`IDAGetSensStats` returns all of the above
   sensitivity-related solver  statistics as a group.

   **Arguments:**
     * ``ida_mem`` -- pointer to the IDAS memory block.
     * ``nresSevals`` -- number of calls to the sensitivity residual function.
     * ``nresevalsS`` -- number of calls to the user-supplied DAE residual
       function for sensitivity evaluations.
     * ``nSetfails`` -- number of error test failures.
     * ``nlinsetupsS`` -- number of calls to the linear solver setup function.

   **Return value:**
     * ``IDA_SUCCESS`` -- The optional output values have been successfully set.
     * ``IDA_MEM_NULL`` -- The ``ida_mem`` pointer is ``NULL``.
     * ``IDA_NO_SENS`` -- Forward sensitivity analysis was not initialized.


.. c:function:: int IDAGetSensErrWeights(void * ida_mem, N_Vector * eSweight)

   The function :c:func:`IDAGetSensErrWeights` returns the sensitivity error weight
   vectors at the current time. These are the reciprocals of the :math:`W_i` of
   :eq:`IDAS_errwt` for the sensitivity variables.

   **Arguments:**
     * ``ida_mem`` -- pointer to the IDAS memory block.
     * ``eSweight`` -- pointer to the array of error weight vectors.

   **Return value:**
     * ``IDA_SUCCESS`` -- The optional output value has been successfully set.
     * ``IDA_MEM_NULL`` -- The ``ida_mem`` pointer is ``NULL``.
     * ``IDA_NO_SENS`` -- Forward sensitivity analysis was not initialized.

   **Notes:**
      The user must allocate memory for ``eweightS``.


.. c:function:: int IDAGetSensNumNonlinSolvIters(void * ida_mem, long int * nSniters)

   The function :c:func:`IDAGetSensNumNonlinSolvIters` returns the  number of
   nonlinear iterations performed for  sensitivity calculations.

   **Arguments:**
     * ``ida_mem`` -- pointer to the IDAS memory block.
     * ``nSniters`` -- number of nonlinear iterations performed.

   **Return value:**
     * ``IDA_SUCCESS`` -- The optional output value has been successfully set.
     * ``IDA_MEM_NULL`` -- The ``ida_mem`` pointer is ``NULL``.
     * ``IDA_NO_SENS`` -- Forward sensitivity analysis was not initialized.
     * ``IDA_MEM_FAIL`` -- The SUNNONLINSOL module is ``NULL``.

   **Notes:**
      This counter is incremented only if ``ism`` was ``IDA_STAGGERED`` or
      in the call to :c:func:`IDASensInit`.


.. c:function:: int IDAGetSensNumNonlinSolvConvFails(void * ida_mem, long int * nSncfails)

   The function :c:func:`IDAGetSensNumNonlinSolvConvFails` returns the  number of
   nonlinear convergence failures that have occurred for  sensitivity
   calculations.

   **Arguments:**
     * ``ida_mem`` -- pointer to the IDAS memory block.
     * ``nSncfails`` -- number of nonlinear convergence failures.

   **Return value:**
     * ``IDA_SUCCESS`` -- The optional output value has been successfully set.
     * ``IDA_MEM_NULL`` -- The ``ida_mem`` pointer is ``NULL``.
     * ``IDA_NO_SENS`` -- Forward sensitivity analysis was not initialized.

   **Notes:**
      This counter is incremented only if ``ism`` was ``IDA_STAGGERED`` or
      in the call to :c:func:`IDASensInit`.


.. c:function:: int IDAGetSensNonlinSolvStats(void * ida_mem, long int * nSniters, long int * nSncfails)

   The function :c:func:`IDAGetSensNonlinSolvStats` returns the sensitivity-related
   nonlinear solver statistics as a group.

   **Arguments:**
     * ``ida_mem`` -- pointer to the IDAS memory block.
     * ``nSniters`` -- number of nonlinear iterations performed.
     * ``nSncfails`` -- number of nonlinear convergence failures.

   **Return value:**
     * ``IDA_SUCCESS`` -- The optional output values have been successfully set.
     * ``IDA_MEM_NULL`` -- The ``ida_mem`` pointer is ``NULL``.
     * ``IDA_NO_SENS`` -- Forward sensitivity analysis was not initialized.
     * ``IDA_MEM_FAIL`` -- The SUNNONLINSOL module is ``NULL``.


.. _IDAS.Usage.FSA.user_callable.optional_output_ic:

Initial condition calculation optional output functions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The sensitivity consistent initial conditions found by IDAS (after a successful
call to :c:func:`IDACalcIC`) can be obtained by calling the following function:

.. c:function:: int IDAGetSensConsistentIC(void * ida_mem, N_Vector * yyS0_mod, N_Vector * ypS0_mod)

   The function :c:func:`IDAGetSensConsistentIC` returns the corrected initial
   conditions calculated by :c:func:`IDACalcIC` for sensitivities variables.

   **Arguments:**
     * ``ida_mem`` -- pointer to the IDAS memory block.
     * ``yyS0_mod`` -- a pointer to an array of ``Ns`` vectors containing consistent sensitivity vectors.
     * ``ypS0_mod`` -- a pointer to an array of ``Ns`` vectors containing consistent sensitivity derivative vectors.

   **Return value:**
     * ``IDA_SUCCESS`` -- :c:func:`IDAGetSensConsistentIC` succeeded.
     * ``IDA_MEM_NULL`` -- The ``ida_mem`` pointer is ``NULL``.
     * ``IDA_NO_SENS`` -- The function :c:func:`IDASensInit` has not been previously called.
     * ``IDA_ILL_INPUT`` -- :c:func:`IDASolve` has been already called.

   **Notes:**
      If the consistent sensitivity vectors or consistent derivative vectors  are not desired, pass ``NULL`` for the corresponding argument.

      .. warning::
         The user must allocate space for ``yyS0_mod`` and ``ypS0_mod``  (if not ``NULL``).


.. _IDAS.Usage.FSA.user_supplied:

User-supplied routines for forward sensitivity analysis
-------------------------------------------------------

In addition to the required and optional user-supplied routines described in
:numref:`IDAS.Usage.SIM.user_supplied`, when using IDAS for forward sensitivity
analysis, the user has the option of providing a routine that calculates the
residual of the sensitivity equations :eq:`IDAS_sens_eqns`.

By default, IDAS uses difference quotient approximation routines for the
residual of the sensitivity equations. However, IDAS allows the option for
user-defined sensitivity residual routines (which also provides a mechanism for
interfacing IDAS to routines generated by automatic differentiation).

The user may provide the residuals of the sensitivity equations :eq:`IDAS_sens_eqns`
for all sensitivity parameters at once, through a function of type
:c:type:`IDASensResFn` defined by:


.. c:type:: int (*IDASensResFn)(int Ns, realtype t, N_Vector yy, N_Vector yp, N_Vector resval, N_Vector *yS, N_Vector *ypS, N_Vector *resvalS, void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)

   This function computes the sensitivity residual for all sensitivity
   equations.  It must compute the vectors
   :math:`\left({\partial F}/{\partial y_i}\right)s_i(t) + \left({\partial F}/{\partial \dot y}\right) \dot{s}_i(t) + \left({\partial F}/{\partial p_i}\right)`
   and store them in ``resvalS[i]``.

   **Arguments:**
     * ``Ns`` -- is the number of sensitivities.
     * ``t`` -- is the current value of the independent variable.
     * ``yy`` -- is the current value of the state vector, :math:`y(t)` .
     * ``yp`` -- is the current value of :math:`\dot{y}(t)` .
     * ``resval`` -- contains the current value :math:`F` of the original DAE residual.
     * ``yS`` -- contains the current values of the sensitivities :math:`s_i` .
     * ``ypS`` -- contains the current values of the sensitivity derivatives :math:`\dot{s}_i` .
     * ``resvalS`` -- contains the output sensitivity residual vectors. Memory allocation for ``resvalS`` is handled within IDAS.
     * ``user_data`` -- is a pointer to user data.
     * ``tmp1``, ``tmp2``, ``tmp3`` -- are ``N_Vector`` s of length :math:`N` which can be used as temporary storage.

   **Return value:**
      An :c:func:`IDASensResFn` should return 0 if successful, a positive value if a recoverable
      error occurred (in which case IDAS will attempt to correct), or a negative
      value if it failed unrecoverably (in which case the integration is halted and
      ``IDA_SRES_FAIL`` is returned).

   **Notes:**
      There is one situation in which recovery is not possible even if
      :c:func:`IDASensResFn` function returns a recoverable error flag.  That is  when
      this occurs at the very first call to the :c:func:`IDASensResFn`, in  which case
      IDAS returns ``IDA_FIRST_RES_FAIL``.


.. _IDAS.Usage.FSA.quad:

Integration of quadrature equations depending on forward sensitivities
----------------------------------------------------------------------

IDAS provides support for integration of quadrature equations that depends not
only on the state variables but also on forward sensitivities.

The following is an overview of the sequence of calls in a user's main program
in this situation. Steps that are changed from the skeleton program presented
in :numref:`IDAS.Usage.SIM.skeleton_sim` are bolded. See also
:numref:`IDAS.Usage.Purequad`.

#. Initialize parallel or multi-threaded environment, if appropriate

#. Create the SUNDIALS context object

#. Set vector of initial values

#. Create matrix object

#. Create linear solver object

#. Set linear solver optional inputs

#. Create nonlinear solver object

#. Create IDAS object

#. Initialize IDAS solver

#. Specify integration tolerances

#. Attach linear solver

#. Set linear solver optional inputs

#. Attach nonlinear solver

#. Set nonlinear solver optional inputs

#. Set sensitivity initial values

#. Activate sensitivity calculations

#. Set sensitivity integration tolerances

#. Create sensitivity nonlinear solver

#. Attach the sensitivity nonlinear solver

#. Set sensitivity nonlinear solver optional inputs

#. **Set vector of initial values for quadrature variables**

   Typically, the quadrature variables should be initialized to :math:`0`.

#. **Initialize sensitivity-dependent quadrature integration**

   Call :c:func:`IDAQuadSensInit` to specify the quadrature equation
   right-hand side function and to allocate internal memory related to
   quadrature integration.

#. Specify rootfinding problem

#. **Set optional inputs**

   Call :c:func:`IDASetQuadSensErrCon` to indicate whether or not quadrature
   variables should be used in the step size control mechanism. If so, one of
   the ``IDAQuadSens*tolerances`` functions must be called to specify the
   integration tolerances for quadrature variables. See
   :numref:`IDAS.Usage.Purequad.quad_optional_input` for details.

#. Correct initial values

#. Advance solution in time

#. Extract sensitivity solution

#. **Extract sensitivity-dependent quadrature variables**

   Call :c:func:`IDAGetQuadSens`, :c:func:`IDAGetQuadSens1`,
   :c:func:`IDAGetQuadSensDky` or :c:func:`IDAGetQuadSensDky1` to obtain the
   values of the quadrature variables or their derivatives at the current time.

#. **Get optional outputs**

   Call ``IDAGetQuadSens*`` functions to obtain optional output related to the
   integration of sensitivity-dependent quadratures. See
   :numref:`IDAS.Usage.FSA.quad.quad_sens_optional_output` for details.

#. Deallocate memory

#. Finalize MPI, if used


.. _IDAS.Usage.FSA.quad.quad_init:

Sensitivity-dependent quadrature initialization and deallocation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The function :c:func:`IDAQuadSensInit` activates integration of quadrature equations
depending on sensitivities and allocates internal memory related to these
calculations. If ``rhsQS`` is input as ``NULL``, then IDAS uses an internal
function that computes difference quotient approximations to the functions
:math:`\bar q_i = (\partial q / \partial y) s_i + (\partial q / \partial \dot{y}) \dot{s}_i + \partial q / \partial p_i`,
in the notation of :eq:`IDAS_QUAD`. The form of the call to this function is as follows:

.. c:function:: int IDAQuadSensInit(void * ida_mem, IDAQuadSensRhsFn rhsQS, N_Vector * yQS0)

   The function :c:func:`IDAQuadSensInit` provides required problem specifications,  allocates internal memory, and initializes quadrature integration.

   **Arguments:**
     * ``ida_mem`` -- pointer to the IDAS memory block returned by :c:func:`IDACreate`.
     * ``rhsQS`` -- is the :c:type:`IDAQuadSensRhsFn` function which computes :math:`f_{QS}` , the right-hand side of the sensitivity-dependent quadrature equations.
     * ``yQS0`` -- contains the initial values of sensitivity-dependent quadratures.

   **Return value:**
     * ``IDA_SUCCESS`` -- The call to :c:func:`IDAQuadSensInit` was successful.
     * ``IDA_MEM_NULL`` -- The IDAS memory was not initialized by a prior call to :c:func:`IDACreate`.
     * ``IDA_MEM_FAIL`` -- A memory allocation request failed.
     * ``IDA_NO_SENS`` -- The sensitivities were not initialized by a prior call to :c:func:`IDASensInit`.
     * ``IDA_ILL_INPUT`` -- The parameter ``yQS0`` is ``NULL``.

   **Notes:**
      .. warning::

         Before calling :c:func:`IDAQuadSensInit`, the user must enable the
         sensitivites  by calling  :c:func:`IDASensInit`.  If an error occurred,
         :c:func:`IDAQuadSensInit` also sends an error message to the  error handler
         function.

In terms of the number of quadrature variables :math:`N_q` and maximum method
order ``maxord``, the size of the real workspace is increased as follows:

* Base value:
  :math:`\text{\texttt{lenrw}} = \text{\texttt{lenrw}} + (\text{\texttt{maxord}} + 5) N_q`

* If :c:func:`IDAQuadSensSVtolerances` is called:
  :math:`\text{\texttt{lenrw}} = \text{\texttt{lenrw}} + N_q N_s`

and the size of the integer workspace is increased as follows:

* Base value:
  :math:`\text{\texttt{leniw}} = \text{\texttt{leniw}} + (\text{\texttt{maxord}} + 5) N_q`

* If :c:func:`IDAQuadSensSVtolerances` is called:
  :math:`\text{\texttt{leniw}} = \text{\texttt{leniw}} + N_q N_s`

The function :c:func:`IDAQuadSensReInit`, useful during the solution of a sequence of
problems of same size, reinitializes the quadrature related internal memory and
must follow a call to :c:func:`IDAQuadSensInit`. The number ``Nq`` of quadratures as
well as the number ``Ns`` of sensitivities are assumed to be unchanged from the
prior call to :c:func:`IDAQuadSensInit`. The call to the :c:func:`IDAQuadSensReInit`
function has the form:

.. c:function:: int IDAQuadSensReInit(void * ida_mem, N_Vector * yQS0)

   The function :c:func:`IDAQuadSensReInit` provides required problem specifications  and reinitializes the sensitivity-dependent quadrature integration.

   **Arguments:**
     * ``ida_mem`` -- pointer to the IDAS memory block.
     * ``yQS0`` -- contains the initial values of sensitivity-dependent quadratures.

   **Return value:**
     * ``IDA_SUCCESS`` -- The call to :c:func:`IDAQuadSensReInit` was successful.
     * ``IDA_MEM_NULL`` -- The IDAS memory was not initialized by a prior call to :c:func:`IDACreate`.
     * ``IDA_NO_SENS`` -- Memory space for the sensitivity calculation was not allocated by a prior call to :c:func:`IDASensInit`.
     * ``IDA_NO_QUADSENS`` -- Memory space for the sensitivity quadratures integration was not allocated by a prior call to :c:func:`IDAQuadSensInit`.
     * ``IDA_ILL_INPUT`` -- The parameter ``yQS0`` is ``NULL``.

   **Notes:**
      If an error occurred, :c:func:`IDAQuadSensReInit` also sends an error message to the  error handler function.


.. c:function:: void IDAQuadSensFree(void* ida_mem);

   The function :c:func:`IDAQuadSensFree` frees the memory allocated for
   sensitivity quadrature integration.

   **Arguments:**
     * ``ida_mem`` -- pointer to the IDAS memory block.

   **Return value:**
      There is no return value.

   **Notes:**
      In general, :c:func:`IDAQuadSensFree` need not be called by the user as it
      is called automatically by :c:func:`IDAFree`.


IDAS solver function
^^^^^^^^^^^^^^^^^^^^

Even if quadrature integration was enabled, the call to the main solver function
:c:func:`IDASolve` is exactly the same as in :numref:`IDAS.Usage.SIM`.
However, in this case the return value ``flag`` can also be one of the
following:

- ``IDA_QSRHS_FAIL`` -- the sensitivity quadrature right-hand side function failed in an unrecoverable manner.

- ``IDA_FIRST_QSRHS_ERR`` -- the sensitivity quadrature right-hand side function failed at the first call.

- ``IDA_REP_QSRHS_ERR`` -- convergence test failures occurred too many times due to repeated
  recoverable errors in the quadrature right-hand side function. The
  ``IDA_REP_RES_ERR`` will also be returned if the quadrature right-hand side
  function had repeated recoverable errors during the estimation of an initial
  step size (assuming the sensitivity quadrature variables are included in the
  error tests).

.. _IDAS.Usage.FSA.quad.quad_sens_get:

Sensitivity-dependent quadrature extraction functions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If sensitivity-dependent quadratures have been initialized by a call to :c:func:`IDAQuadSensInit`, or reinitialized by a call
to :c:func:`IDAQuadSensReInit`, then IDAS computes a solution, sensitivities, and quadratures depending on sensitivities
at time ``t``. However, :c:func:`IDASolve` will still return only the solutions :math:`y` and :math:`\dot{y}`.
Sensitivity-dependent quadratures can be obtained using one of the following
functions:

.. c:function:: int IDAGetQuadSens(void * ida_mem, realtype * tret, N_Vector * yQS)

   The function :c:func:`IDAGetQuadSens` returns the quadrature sensitivity  solution vectors after a successful return from :c:func:`IDASolve`.

   **Arguments:**
     * ``ida_mem`` -- pointer to the memory previously allocated by :c:func:`IDAInit`.
     * ``tret`` -- the time reached by the solver output.
     * ``yQS`` -- array of ``Ns`` computed sensitivity-dependent quadrature vectors. This array of vectors must be allocated by the user.

   **Return value:**
     * ``IDA_SUCCESS`` -- :c:func:`IDAGetQuadSens` was successful.
     * ``IDA_MEM_NULL`` -- ``ida_mem`` was NULL.
     * ``IDA_NO_SENS`` -- Sensitivities were not activated.
     * ``IDA_NO_QUADSENS`` -- Quadratures depending on the sensitivities were not activated.
     * ``IDA_BAD_DKY`` -- ``yQS`` or one of the ``yQS[i]`` is ``NULL``.


The function :c:func:`IDAGetQuadSensDky` computes the ``k``-th derivatives of
the interpolating polynomials for the sensitivity-dependent quadrature variables
at time ``t``. This function is called by :c:func:`IDAGetQuadSens` with
``k = 0``, but may also be called directly by the user.

.. c:function:: int IDAGetQuadSensDky(void * ida_mem, realtype t, int k, N_Vector* dkyQS)

   The function :c:func:`IDAGetQuadSensDky` returns derivatives of the quadrature sensitivities  solution vectors after a successful return from :c:func:`IDASolve`.

   **Arguments:**
     * ``ida_mem`` -- pointer to the memory previously allocated by :c:func:`IDAInit`.
     * ``t`` -- the time at which information is requested. The time ``t`` must fall within the interval defined by the last successful step taken by IDAS.
     * ``k`` -- order of the requested derivative. ``k`` must be in the range :math:`0, 1, ..., klast` where :math:`klast` is the method order of the last successful step.
     * ``dkyQS`` -- array of ``Ns`` vectors containing the derivatives. This vector array must be allocated by the user.

   **Return value:**
     * ``IDA_SUCCESS`` -- :c:func:`IDAGetQuadSensDky` succeeded.
     * ``IDA_MEM_NULL`` -- ``ida_mem`` was ``NULL``.
     * ``IDA_NO_SENS`` -- Sensitivities were not activated.
     * ``IDA_NO_QUADSENS`` -- Quadratures depending on the sensitivities were not activated.
     * ``IDA_BAD_DKY`` -- ``dkyQS`` or one of the vectors ``dkyQS[i]`` is ``NULL``.
     * ``IDA_BAD_K`` -- ``k`` is not in the range :math:`0, 1, ..., klast`.
     * ``IDA_BAD_T`` -- The time ``t`` is not in the allowed range.


Quadrature sensitivity solution vectors can also be extracted separately for
each parameter in turn through the functions ``IDAGetQuadSens1`` and
``IDAGetQuadSensDky1``, defined as follows:

.. c:function:: int IDAGetQuadSens1(void * ida_mem, realtype * tret, int is, N_Vector yQS)

   The function ``IDAGetQuadSens1`` returns the ``is``-th sensitivity  of quadratures after a successful return from :c:func:`IDASolve`.

   **Arguments:**
     * ``ida_mem`` -- pointer to the memory previously allocated by :c:func:`IDAInit`.
     * ``tret`` -- the time reached by the solver output.
     * ``is`` -- specifies which sensitivity vector is to be returned :math:`0\le` ``is`` :math:`< N_s`.
     * ``yQS`` -- the computed sensitivity-dependent quadrature vector. This vector must be allocated by the user.

   **Return value:**
     * ``IDA_SUCCESS`` -- ``IDAGetQuadSens1`` was successful.
     * ``IDA_MEM_NULL`` -- ``ida_mem`` was ``NULL``.
     * ``IDA_NO_SENS`` -- Forward sensitivity analysis was not initialized.
     * ``IDA_NO_QUADSENS`` -- Quadratures depending on the sensitivities were not activated.
     * ``IDA_BAD_IS`` -- The index ``is`` is not in the allowed range.
     * ``IDA_BAD_DKY`` -- ``yQS`` is ``NULL``.


.. c:function:: int IDAGetQuadSensDky1(void * ida_mem, realtype t, int k, int is, N_Vector dkyQS)

   The function ``IDAGetQuadSensDky1`` returns the ``k``-th derivative of the  ``is``-th sensitivity solution vector after a successful  return from :c:func:`IDASolve`.

   **Arguments:**
     * ``ida_mem`` -- pointer to the memory previously allocated by :c:func:`IDAInit`.
     * ``t`` -- specifies the time at which sensitivity information is requested. The time ``t`` must fall within the interval defined by the last successful step taken by IDAS.
     * ``k`` -- order of derivative. ``k`` must be in the range :math:`0, 1, ..., klast` where :math:`klast` is the method order of the last successful step.
     * ``is`` -- specifies the sensitivity derivative vector to be returned :math:`0\le` ``is`` :math:`< N_s`.
     * ``dkyQS`` -- the vector containing the derivative. The space for ``dkyQS`` must be allocated by the user.

   **Return value:**
     * ``IDA_SUCCESS`` -- ``IDAGetQuadDky1`` succeeded.
     * ``IDA_MEM_NULL`` -- ``ida_mem`` was ``NULL``.
     * ``IDA_NO_SENS`` -- Forward sensitivity analysis was not initialized.
     * ``IDA_NO_QUADSENS`` -- Quadratures depending on the sensitivities were not activated.
     * ``IDA_BAD_DKY`` -- ``dkyQS`` is ``NULL``.
     * ``IDA_BAD_IS`` -- The index ``is`` is not in the allowed range.
     * ``IDA_BAD_K`` -- ``k`` is not in the range :math:`0, 1, ..., klast`.
     * ``IDA_BAD_T`` -- The time ``t`` is not in the allowed range.


.. _IDAS.Usage.FSA.quad.quad_sens_optional_input:

Optional inputs for sensitivity-dependent quadrature integration
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

IDAS provides the following optional input functions to control the integration
of sensitivity-dependent quadrature equations.

.. c:function:: int IDASetQuadSensErrCon(void * ida_mem, booleantype errconQS)

   The function :c:func:`IDASetQuadSensErrCon` specifies whether or not the  quadrature variables are to be used in the local error control mechanism.  If they are, the user must specify the error tolerances for the quadrature  variables by calling :c:func:`IDAQuadSensSStolerances`,  :c:func:`IDAQuadSensSVtolerances`, or :c:func:`IDAQuadSensEEtolerances`.

   **Arguments:**
     * ``ida_mem`` -- pointer to the IDAS memory block.
     * ``errconQS`` -- specifies whether sensitivity quadrature variables are included ``SUNTRUE`` or not ``SUNFALSE`` in the error control mechanism.

   **Return value:**
     * ``IDA_SUCCESS`` -- The optional value has been successfully set.
     * ``IDA_MEM_NULL`` -- The ``ida_mem`` pointer is ``NULL``.
     * ``IDA_NO_SENS`` -- Sensitivities were not activated.
     * ``IDA_NO_QUADSENS`` -- Quadratures depending on the sensitivities were not activated.

   **Notes:**
      By default, ``errconQS`` is set to ``SUNFALSE``.

      .. warning::
         It is illegal to call :c:func:`IDASetQuadSensErrCon` before a call  to :c:func:`IDAQuadSensInit`.


If the quadrature variables are part of the step size control mechanism, one of
the following functions must be called to specify the integration tolerances for
quadrature variables.

.. c:function:: int IDAQuadSensSStolerances(void * ida_mem, realtype reltolQS, realtype* abstolQS)

   The function :c:func:`IDAQuadSensSStolerances` specifies scalar relative and absolute  tolerances.

   **Arguments:**
     * ``ida_mem`` -- pointer to the IDAS memory block.
     * ``reltolQS`` --  tolerances is the scalar relative error tolerance.
     * ``abstolQS`` -- is a pointer to an array containing the ``Ns`` scalar absolute error tolerances.

   **Return value:**
     * ``IDA_SUCCESS`` -- The optional value has been successfully set.
     * ``IDA_MEM_NULL`` -- The ``ida_mem`` pointer is ``NULL``.
     * ``IDA_NO_SENS`` -- Sensitivities were not activated.
     * ``IDA_NO_QUADSENS`` -- Quadratures depending on the sensitivities were not activated.
     * ``IDA_ILL_INPUT`` -- One of the input tolerances was negative.


.. c:function:: int IDAQuadSensSVtolerances(void * ida_mem, realtype reltolQS, N_Vector* abstolQS)

   The function :c:func:`IDAQuadSensSVtolerances` specifies scalar relative and  vector absolute tolerances.

   **Arguments:**
     * ``ida_mem`` -- pointer to the IDAS memory block.
     * ``reltolQS`` --  tolerances is the scalar relative error tolerance.
     * ``abstolQS`` -- is an array of ``Ns`` variables of type ``N_Vector``. The ``N_Vector`` from ``abstolS[is]`` specifies the vector tolerances for ``is`` -th quadrature sensitivity.

   **Return value:**
     * ``IDA_SUCCESS`` -- The optional value has been successfully set.
     * ``IDA_NO_QUAD`` -- Quadrature integration was not initialized.
     * ``IDA_MEM_NULL`` -- The ``ida_mem`` pointer is ``NULL``.
     * ``IDA_NO_SENS`` -- Sensitivities were not activated.
     * ``IDA_NO_QUADSENS`` -- Quadratures depending on the sensitivities were not activated.
     * ``IDA_ILL_INPUT`` -- One of the input tolerances was negative.


.. c:function:: int IDAQuadSensEEtolerances(void * ida_mem)

   The function :c:func:`IDAQuadSensEEtolerances` specifies that the tolerances
   for  the sensitivity-dependent quadratures should be estimated from those
   provided  for the pure quadrature variables.

   **Arguments:**
     * ``ida_mem`` -- pointer to the IDAS memory block.

   **Return value:**
     * ``IDA_SUCCESS`` -- The optional value has been successfully set.
     * ``IDA_MEM_NULL`` -- The ``ida_mem`` pointer is ``NULL``.
     * IDA_NO_SENS -- Sensitivities were not activated.
     * ``IDA_NO_QUADSENS`` -- Quadratures depending on the sensitivities were not activated.

   **Notes:**
      When :c:func:`IDAQuadSensEEtolerances`  is used, before calling
      :c:func:`IDASolve`,  integration of pure quadratures must be initialized
      (see :numref:`IDAS.Usage.Purequad`)  and tolerances for pure
      quadratures must be also specified  (see
      :numref:`IDAS.Usage.Purequad.quad_optional_input`).


.. _IDAS.Usage.FSA.quad.quad_sens_optional_output:

Optional outputs for sensitivity-dependent quadrature integration
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

IDAS provides the following functions that can be used to obtain solver
performance information related to quadrature integration.

.. c:function:: int IDAGetQuadSensNumRhsEvals(void * ida_mem, long int * nrhsQSevals)

   The function :c:func:`IDAGetQuadSensNumRhsEvals` returns the  number of calls
   made to the user's quadrature right-hand side function.

   **Arguments:**
     * ``ida_mem`` -- pointer to the IDAS memory block.
     * ``nrhsQSevals`` -- number of calls made to the user's ``rhsQS`` function.

   **Return value:**
     * ``IDA_SUCCESS`` -- The optional output value has been successfully set.
     * ``IDA_MEM_NULL`` -- The ``ida_mem`` pointer is ``NULL``.
     * ``IDA_NO_QUADSENS`` -- Sensitivity-dependent quadrature integration has not been initialized.


.. c:function:: int IDAGetQuadSensNumErrTestFails(void * ida_mem, long int * nQSetfails)

   The function :c:func:`IDAGetQuadSensNumErrTestFails` returns the  number of
   local error test failures due to quadrature variables.

   **Arguments:**
     * ``ida_mem`` -- pointer to the IDAS memory block.
     * ``nQSetfails`` -- number of error test failures due to quadrature variables.

   **Return value:**
     * ``IDA_SUCCESS`` -- The optional output value has been successfully set.
     * ``IDA_MEM_NULL`` -- The ``ida_mem`` pointer is ``NULL``.
     * ``IDA_NO_QUADSENS`` -- Sensitivity-dependent quadrature integration has not been initialized.


.. c:function:: int IDAGetQuadSensErrWeights(void * ida_mem, N_Vector * eQSweight)

   The function :c:func:`IDAGetQuadSensErrWeights` returns the quadrature error
   weights  at the current time.

   **Arguments:**
     * ``ida_mem`` -- pointer to the IDAS memory block.
     * ``eQSweight`` -- array of quadrature error weight vectors at the current time.

   **Return value:**
     * ``IDA_SUCCESS`` -- The optional output value has been successfully set.
     * ``IDA_MEM_NULL`` -- The ``ida_mem`` pointer is ``NULL``.
     * ``IDA_NO_QUADSENS`` -- Sensitivity-dependent quadrature integration has not been initialized.

   **Notes:**
      .. warning::
         The user must allocate memory for ``eQSweight``.  If quadratures were
         not included in the error control mechanism (through a  call to
         :c:func:`IDASetQuadSensErrCon` with ``errconQS=SUNTRUE``),
         :c:func:`IDAGetQuadSensErrWeights` does not set the ``eQSweight``
         vector.


.. c:function:: int IDAGetQuadSensStats(void * ida_mem, long int * nrhsQSevals, long int * nQSetfails)

   The function :c:func:`IDAGetQuadSensStats` returns the IDAS integrator
   statistics  as a group.

   **Arguments:**
     * ``ida_mem`` -- pointer to the IDAS memory block.
     * ``nrhsQSevals`` -- number of calls to the user's ``rhsQS`` function.
     * ``nQSetfails`` -- number of error test failures due to quadrature variables.

   **Return value:**
     * ``IDA_SUCCESS`` -- the optional output values have been successfully set.
     * ``IDA_MEM_NULL`` -- the ``ida_mem`` pointer is ``NULL``.
     * ``IDA_NO_QUADSENS`` -- Sensitivity-dependent quadrature integration has not been initialized.


.. _IDAS.Usage.FSA.quad.user_supplied:

User-supplied function for sensitivity-dependent quadrature integration
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

For the integration of sensitivity-dependent quadrature equations, the user must
provide a function that defines the residual of those quadrature equations. For
the sensitivities of quadratures :eq:`IDAS_QUAD` with integrand :math:`q`, the
appropriate residual functions are given by :math:`\bar{q}_i = {\partial
q}/{\partial y} s_i + {\partial q}/{\partial \dot{y}} \dot{s}_i + {\partial
q}{\partial p_i}`.  This user function must be of type
:c:type:`IDAQuadSensRhsFn` defined as follows:

.. c:type:: int (*IDAQuadSensRhsFn)(int Ns, realtype t, N_Vector yy, N_Vector yp, N_Vector *yyS, N_Vector *ypS, N_Vector rrQ, N_Vector *rhsvalQS, void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)

   This function computes the sensitivity quadrature equation right-hand side
   for a given value  of the independent variable :math:`t` and state vector
   :math:`y`.

   **Arguments:**
     * ``Ns`` -- is the number of sensitivity vectors.
     * ``t`` -- is the current value of the independent variable.
     * ``yy`` -- is the current value of the dependent variable vector, :math:`y(t)`.
     * ``yp`` -- is the current value of the dependent variable vector, :math:`\dot{y}(t)`.
     * ``yyS`` -- is an array of ``Ns`` variables of type ``N_Vector``
       containing the dependent sensitivity vectors :math:`s_i`.
     * ``ypS`` -- is an array of ``Ns`` variables of type ``N_Vector`` containing the dependent sensitivity derivatives :math:`\dot{s}_i`.
     * ``rrQ`` -- is the current value of the quadrature right-hand side :math:`q`.
     * ``rhsvalQS`` -- contains the ``Ns`` output vectors.
     * ``user_data`` -- is the ``user_data`` pointer passed to :c:func:`IDASetUserData`.
     * ``tmp1``, ``tmp2``, ``tmp3`` -- are ``N_Vector`` s which can be used as temporary storage.

   **Return value:**
      An :c:type:`IDAQuadSensRhsFn` should return 0 if successful, a positive
      value if a recoverable error occurred (in which case IDAS will attempt to
      correct), or a negative value if it failed unrecoverably (in which case
      the integration is halted and ``IDA_QRHS_FAIL`` is returned).

   **Notes:**

   Allocation of memory for ``rhsvalQS`` is automatically handled within IDAS.

   Both ``yy`` and ``yp`` are of type ``N_Vector`` and both ``yyS`` and ``ypS``
   are pointers to an array containing ``Ns`` vectors of type ``N_Vector``.  It
   is the user's responsibility to access the vector data consistently
   (including the use of the correct accessor macros from each ``N_Vector``
   implementation).

   There is one situation in which recovery is not possible even if
   :c:type:`IDAQuadSensRhsFn` function returns a recoverable error flag.  That
   is when this occurs at the very first call to the :c:type:`IDAQuadSensRhsFn`,
   in which case IDAS returns ``IDA_FIRST_QSRHS_ERR``).


.. _IDAS.Usage.FSA.partial:

Note on using partial error control
-----------------------------------

For some problems, when sensitivities are excluded from the error control test,
the behavior of IDAS may appear at first glance to be erroneous. One would
expect that, in such cases, the sensitivity variables would not influence in any
way the step size selection.

The short explanation of this behavior is that the step size selection
implemented by the error control mechanism in IDAS is based on the magnitude of
the correction calculated by the nonlinear solver. As mentioned in
:numref:`IDAS.Usage.FSA.user_callable.sensi_init`, even with partial error
control selected in the call to :c:func:`IDASensInit`, the sensitivity variables are
included in the convergence tests of the nonlinear solver.

When using the simultaneous corrector method :numref:`IDAS.Mathematics.FSA`, the
nonlinear system that is solved at each step involves both the state and
sensitivity equations. In this case, it is easy to see how the sensitivity
variables may affect the convergence rate of the nonlinear solver and therefore
the step size selection. The case of the staggered corrector approach is more
subtle. The sensitivity variables at a given step are computed only once the
solver for the nonlinear state equations has converged. However, if the
nonlinear system corresponding to the sensitivity equations has convergence
problems, IDAS will attempt to improve the initial guess by reducing the step
size in order to provide a better prediction of the sensitivity variables.
Moreover, even if there are no convergence failures in the solution of the
sensitivity system, IDAS may trigger a call to the linear solver's setup routine
which typically involves reevaluation of Jacobian information (Jacobian
approximation in the case of matrix-based linear solvers, or preconditioner data in the
case of the Krylov solvers). The new Jacobian information will be used by
subsequent calls to the nonlinear solver for the state equations and, in this
way, potentially affect the step size selection.

When using the simultaneous corrector method it is not possible to decide
whether nonlinear solver convergence failures or calls to the linear solver
setup routine have been triggered by convergence problems due to the state or
the sensitivity equations. When using one of the staggered corrector methods,
however, these situations can be identified by carefully monitoring the
diagnostic information provided through optional outputs. If there are no
convergence failures in the sensitivity nonlinear solver, and none of the calls
to the linear solver setup routine were made by the sensitivity nonlinear
solver, then the step size selection is not affected by the sensitivity
variables.

Finally, the user must be warned that the effect of appending sensitivity
equations to a given system of DAEs on the step size selection (through the
mechanisms described above) is problem-dependent and can therefore lead to
either an increase or decrease of the total number of steps that IDAS takes to
complete the simulation. At first glance, one would expect that the impact of
the sensitivity variables, if any, would be in the direction of increasing the
step size and therefore reducing the total number of steps. The argument for
this is that the presence of the sensitivity variables in the convergence test
of the nonlinear solver can only lead to additional iterations (and therefore a
smaller iteration error), or to additional calls to the linear solver setup
routine (and therefore more up-to-date Jacobian information), both of which will
lead to larger steps being taken by IDAS. However, this is true only locally.
Overall, a larger integration step taken at a given time may lead to step size
reductions at later times, due to either nonlinear solver convergence failures
or error test failures.
