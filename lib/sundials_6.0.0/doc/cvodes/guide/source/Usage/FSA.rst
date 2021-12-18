.. ----------------------------------------------------------------
   SUNDIALS Copyright Start
   Copyright (c) 2002-2021, Lawrence Livermore National Security
   and Southern Methodist University.
   All rights reserved.

   See the top-level LICENSE and NOTICE files for details.

   SPDX-License-Identifier: BSD-3-Clause
   SUNDIALS Copyright End
   ----------------------------------------------------------------

.. _CVODES.Usage.FSA:

Using CVODES for Forward Sensitivity Analysis
=============================================

This chapter describes the use of CVODES to compute solution sensitivities using
forward sensitivity analysis. One of our main guiding principles was to design
the CVODES user interface for forward sensitivity analysis as an extension of
that for IVP integration. Assuming a user main program and user-defined support
routines for IVP integration have already been defined, in order to perform
forward sensitivity analysis the user only has to insert a few more calls into
the main program and (optionally) define an additional routine which computes
the right-hand side of the sensitivity systems :eq:`CVODES_sens_eqns`. The only
departure from this philosophy is due to the :c:type:`CVRhsFn` type definition.
Without changing the definition of this type, the only way to pass values of the
problem parameters to the ODE right-hand side function is to require the user
data structure ``f_data`` to contain a pointer to the array of real parameters
:math:`p`.

CVODES uses various constants for both input and output. These are defined as
needed in this chapter, but for convenience are also listed separately in
:numref:`CVODES.Constants`.

We begin with a brief overview, in the form of a skeleton user program.
Following that are detailed descriptions of the interface to the various
user-callable routines and of the user-supplied routines that were not already
described in :numref:`CVODES.Usage.SIM` or :numref:`CVODES.Usage.purequad`.

.. _CVODES.Usage.FSA.skeleton_sim:

A skeleton of the user’s main program
-------------------------------------

The following is a skeleton of the user’s main program (or calling program) as
an application of CVODES. The user program is to have these steps in the order
indicated, unless otherwise noted. For the sake of brevity, we defer many of the
details to the later sections. As in :numref:`CVODES.Usage.SIM.skeleton_sim`,
most steps are independent of the ``N_Vector``, ``SUNMatrix``,
``SUNLinearSolver``, and ``SUNNonlinearSolver`` implementations used. For the
steps that are not, refer to Chapters :numref:`NVectors`, :numref:`SUNMatrix`,
:numref:`SUNLinSol`, :numref:`SUNNonlinSol` for the specific name of the
function to be called or macro to be referenced.

Differences between the user main program in
:numref:`CVODES.Usage.SIM.skeleton_sim` and the one below start only at step 16.
Steps that are unchanged from the skeleton program presented in
:numref:`CVODES.Usage.SIM.skeleton_sim` are left unbolded.

First, note that no additional header files need be included for forward
sensitivity analysis beyond those for IVP solution
:numref:`CVODES.Usage.SIM.skeleton_sim`.

   #. Initialize parallel or multi-threaded environment, if appropriate

   #. Create the SUNDIALS context object

   #. Set problem dimensions etc.

   #. Set vector of initial values

   #. Create CVODE object

   #. Initialize CVODE solver

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

   #. **Define the sensitivity problem**

      -  **Number of sensitivities** (*required*)
         Set ``Ns`` :math:`= N_s`, the number of parameters with respect to which sensitivities are to be computed.

      -  **Problem parameters** (*optional*)
         If CVODES is to evaluate the right-hand sides of the sensitivity systems, set ``p``, an array of ``Np`` real
         parameters upon which the IVP depends. Only parameters with respect to which sensitivities are (potentially) desired
         need to be included. Attach ``p`` to the user data structure ``user_data``. For example, ``user_data->p = p;``

         If the user provides a function to evaluate the sensitivity right-hand side, ``p`` need not be specified.

      -  **Parameter list** (*optional*)
         If CVODES is to evaluate the right-hand sides of the sensitivity systems, set ``plist``, an array of ``Ns``
         integers to specify the parameters ``p`` with respect to which solution sensitivities are to be computed. If
         sensitivities with respect to the :math:`j`-th parameter ``p[j]`` are desired :math:`(0 \leq` ``j`` :math:`<`
         ``Np``), set :math:`{\text{plist}}_i = j`, for some :math:`i = 0,\ldots,N_s-1`.

         If ``plist`` is not specified, CVODES will compute sensitivities with respect to the first ``Ns`` parameters;
         i.e., :math:`{\text{plist}}_i = i` :math:`(i = 0,\ldots, N_s - 1)`.

         If the user provides a function to evaluate the sensitivity right-hand side, ``plist`` need not be specified.

      -  **Parameter scaling factors** (*optional*)
         If CVODES is to estimate tolerances for the sensitivity solution vectors (based on tolerances for the state
         solution vector) or if CVODES is to evaluate the right-hand sides of the sensitivity systems using the internal
         difference-quotient function, the results will be more accurate if order of magnitude information is provided.

         Set ``pbar``, an array of ``Ns`` positive scaling factors. Typically, if :math:`p_i \ne 0`, the value
         :math:`{\bar p}_i = |p_{\text{plist}_i}|` can be used.

         If ``pbar`` is not specified, CVODES will use :math:`{\bar p}_i = 1.0`.

         If the user provides a function to evaluate the sensitivity right-hand side and specifies tolerances for the
         sensitivity variables, ``pbar`` need not be specified.

         Note that the names for ``p``, ``pbar``, ``plist``, as well as the field *p* of ``user_data`` are arbitrary, but they
         must agree with the arguments passed to :c:func:`CVodeSetSensParams` below.

   #. **Set sensitivity initial conditions**

      Set the ``Ns`` vectors ``yS0[i]`` of initial values for sensitivities (for
      :math:`i=0,\ldots,` ``Ns`` :math:`- 1`), using the appropriate functions
      defined by the particular ``N_Vector`` implementation chosen.

      First, create an array of ``Ns`` vectors by making the appropriate call

      ``yS0 = N_VCloneVectorArray_***(Ns, y0);``

      or

      ``yS0 = N_VCloneVectorArrayEmpty_***(Ns, y0);``

      Here the argument ``y0`` serves only to provide the ``N_Vector`` type for
      cloning.

      Then, for each :math:`i = 0,\ldots,`\ ``Ns`` :math:`- 1`, load initial
      values for the i-th sensitivity vector ``yS0[i]``.

   #. **Activate sensitivity calculations**

      Call :c:func:`CVodeSensInit` or :c:func:`CVodeSensInit1` to activate
      forward sensitivity computations and allocate internal memory for CVODES
      related to sensitivity calculations.

   #. **Set sensitivity tolerances**

      Call :c:func:`CVodeSensSStolerances`, :c:func:`CVodeSensSVtolerances` or
      :c:func:`CVodeEEtolerances`.

   #. **Set sensitivity analysis optional inputs**

      Call ``CVodeSetSens*`` routines to change from their default values any
      optional inputs that control the behavior of CVODES in computing forward
      sensitivities. See :numref:`CVODES.Usage.FSA.user_callable.optional_inputs` for details.

   #. **Create sensitivity nonlinear solver object**

      If using a non-default nonlinear solver (see
      :numref:`CVODES.Usage.FSA.user_callable.nonlin_solv_init`), then create the desired
      nonlinear solver object by calling the appropriate constructor function
      defined by the particular ``SUNNonlinearSolver`` implementation e.g.,

      .. code-block:: c

         NLSSens = SUNNonlinSol_***Sens(...);

      for the ``CV_SIMULTANEOUS`` or ``CV_STAGGERED`` options or

      .. code-block:: c

         NLSSens = SUNNonlinSol_***(...);

      for the ``CV_STAGGERED1`` option where ``***`` is the name of the
      nonlinear solver and ``...`` are constructor specific arguments (see
      :numref:`SUNNonlinSol` for details).

   #. **Attach the sensitivity nonlinear solver module**

      If using a non-default nonlinear solver, then initialize the nonlinear
      solver interface by attaching the nonlinear solver object by calling
      :c:func:`CVodeSetNonlinearSolverSensSim` when using the
      ``CV_SIMULTANEOUS`` corrector method,
      :c:func:`CVodeSetNonlinearSolverSensStg` when using the ``CV_STAGGERED``
      corrector method, or :c:func:`CVodeSetNonlinearSolverSensStg1` when using
      the ``CV_STAGGERED1`` corrector method (see
      :numref:`CVODES.Usage.FSA.user_callable.nonlin_solv_init` for details).


   #. **Set sensitivity nonlinear solver optional inputs**

      Call the appropriate set functions for the selected nonlinear solver
      module to change optional inputs specific to that nonlinear solver. These
      *must* be called after :c:func:`CVodeSensInit` if using the default nonlinear
      solver or after attaching a new nonlinear solver to CVODES, otherwise the
      optional inputs will be overridden by CVODE defaults. See
      :numref:`SUNNonlinSol` for more information on optional inputs.

   #. Specify rootfinding problem (*optional*)

   #. Advance solution in time

   #. **Extract sensitivity solution**

      After each successful return from :c:func:`CVode`, the solution of the
      original IVP is available in the ``y`` argument of :c:func:`CVode`, while
      the sensitivity solution can be extracted into ``yS`` (which can be the
      same as ``yS0``) by calling one of the routines :c:func:`CVodeGetSens`,
      :c:func:`CVodeGetSens1`, :c:func:`CVodeGetSensDky`, or
      :c:func:`CVodeGetSensDky1`.

   #. Get optional outputs

   #. Deallocate memory for solution vector

   #. **Deallocate memory for sensitivity vectors**

      Upon completion of the integration, deallocate memory for the vectors
      ``yS0`` using the appropriate destructor: ``N_VDestroyVectorArray_***(yS0, Ns);``

      If ``yS`` was created from ``realtype`` arrays ``yS_i``, it is the user’s
      responsibility to also free the space for the arrays ``yS0_i``.

   #. Free solver memory

   #. Free nonlinear solver memory (*optional*)

   #. Free linear solver and matrix memory

   #. Free the SUNContext object

   #. Finalize MPI, if used


.. _CVODES.Usage.FSA.user_callable:

User-callable routines for forward sensitivity analysis
-------------------------------------------------------

This section describes the CVODES functions, in addition to those presented in
:numref:`CVODES.Usage.SIM.user_callable`, that are called by the user
to setup and solve a forward sensitivity problem.

.. _CVODES.Usage.FSA.user_callable.sensi_malloc:

Forward sensitivity initialization and deallocation functions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Activation of forward sensitivity computation is done by calling
:c:func:`CVodeSensInit` or :c:func:`CVodeSensInit1`, depending on whether the
sensitivity right-hand side function returns all sensitivities at once or one by
one, respectively. The form of the call to each of these routines is as follows:

.. c:function:: int CVodeSensInit(void * cvode_mem, int Ns, int ism, CVSensRhsFn fS, N_Vector * yS0)

   The routine :c:func:`CVodeSensInit` activates forward sensitivity computations and
   allocates internal memory related to sensitivity calculations.

   **Arguments:**
     * ``cvode_mem`` -- pointer to the CVODES memory block returned by :c:func:`CVodeCreate`.
     * ``Ns`` -- the number of sensitivities to be computed.
     * ``ism`` --  forward sensitivity analysis!correction strategies a flag used to select the sensitivity solution method. Its value can be ``CV_SIMULTANEOUS`` or ``CV_STAGGERED`` :

       * In the ``CV_SIMULTANEOUS`` approach, the state and sensitivity variables are corrected at the same time. If the default Newton nonlinear solver is used, this amounts to performing a modified Newton iteration on the combined nonlinear system;
       * In the ``CV_STAGGERED`` approach, the correction step for the sensitivity variables takes place at the same time for all sensitivity equations, but only after the correction of the state variables has converged and the state variables have passed the local error test;

     * ``fS`` -- is the C function which computes all sensitivity ODE right-hand sides at the same time. For full details see :c:type:`CVSensRhsFn`.
     * ``yS0`` -- a pointer to an array of ``Ns`` vectors containing the initial values of the sensitivities.

   **Return value:**
     * ``CV_SUCCESS`` -- The call to :c:func:`CVodeSensInit` was successful.
     * ``CV_MEM_NULL`` -- The CVODES memory block was not initialized through a previous call to :c:func:`CVodeCreate`.
     * ``CV_MEM_FAIL`` -- A memory allocation request has failed.
     * ``CV_ILL_INPUT`` -- An input argument to :c:func:`CVodeSensInit` has an illegal value.

   **Notes:**
      Passing ``fs == NULL`` indicates using the default internal difference
      quotient sensitivity right-hand side routine.  If an error occurred,
      :c:func:`CVodeSensInit` also sends an error message to the  error handler
      function.

      .. warning::
         It is illegal here to use ``ism = CV_STAGGERED1``. This option
         requires a different type for ``fS`` and can therefore only be used
         with  :c:func:`CVodeSensInit1` (see below).


.. c:function:: int CVodeSensInit1(void * cvode_mem, int Ns, int ism, CVSensRhs1Fn fS1, N_Vector * yS0)

   The routine :c:func:`CVodeSensInit1` activates forward sensitivity computations and
   allocates internal memory related to sensitivity calculations.

   **Arguments:**
     * ``cvode_mem`` -- pointer to the CVODES memory block returned by :c:func:`CVodeCreate`.
     * ``Ns`` -- the number of sensitivities to be computed.
     * ``ism`` --  forward sensitivity analysis!correction strategies a flag used to select the sensitivity solution method. Its value can be ``CV_SIMULTANEOUS`` , ``CV_STAGGERED`` , or ``CV_STAGGERED1`` :

       * In the ``CV_SIMULTANEOUS`` approach, the state and sensitivity variables are corrected at the same time. If the default Newton nonlinear solver is used, this amounts to performing a modified Newton iteration on the combined nonlinear system;
       * In the ``CV_STAGGERED`` approach, the correction step for the sensitivity variables takes place at the same time for all sensitivity equations, but only after the correction of the state variables has converged and the state variables have passed the local error test;
       * In the ``CV_STAGGERED1`` approach, all corrections are done sequentially, first for the state variables and then for the sensitivity variables, one parameter at a time. If the sensitivity variables are not included in the error control, this approach is equivalent to ``CV_STAGGERED``. Note that the ``CV_STAGGERED1`` approach can be used only if the user-provided sensitivity right-hand side function is of type :c:type:`CVSensRhs1Fn`.

     * ``fS1`` -- is the C function which computes the right-hand sides of the sensitivity ODE, one at a time. For full details see :c:type:`CVSensRhs1Fn`.
     * ``yS0`` -- a pointer to an array of ``Ns`` vectors containing the initial values of the sensitivities.

   **Return value:**
     * ``CV_SUCCESS`` -- The call to :c:func:`CVodeSensInit1` was successful.
     * ``CV_MEM_NULL`` -- The CVODES memory block was not initialized through a previous call to :c:func:`CVodeCreate`.
     * ``CV_MEM_FAIL`` -- A memory allocation request has failed.
     * ``CV_ILL_INPUT`` -- An input argument to :c:func:`CVodeSensInit1` has an illegal value.

   **Notes:**
      Passing ``fS1 = NULL`` indicates using the default internal difference
      quotient sensitivity right-hand side routine.  If an error occurred,
      :c:func:`CVodeSensInit1` also sends an error message to the  error handler
      funciton.


In terms of the problem size :math:`N`, number of sensitivity vectors
:math:`N_s`, and maximum method order ``maxord``, the size of the real workspace
is increased as follows:

-  Base value: :math:`\texttt{lenrw} = \texttt{lenrw} + (\texttt{maxord}+5)N_s N`

-  With :c:func:`CVodeSensSVtolerances`: :math:`\texttt{lenrw} = \texttt{lenrw} + N_s N`

the size of the integer workspace is increased as follows:

-  Base value: :math:`\texttt{leniw} = \texttt{leniw} + (\texttt{maxord}+5) N_s N_i`

-  With :c:func:`CVodeSensSVtolerances`: :math:`\texttt{leniw} = \texttt{leniw} + N_s N_i`

where :math:`N_i` is the number of integers in one ``N_Vector``.

The routine :c:func:`CVodeSensReInit`, useful during the solution of a sequence of
problems of same size, reinitializes the sensitivity-related internal memory.
The call to it must follow a call to :c:func:`CVodeSensInit` or :c:func:`CVodeSensInit1`
(and maybe a call to :c:func:`CVodeReInit`). The number ``Ns`` of sensitivities is
assumed to be unchanged since the call to the initialization function. The call
to the :c:func:`CVodeSensReInit` function has the form:

.. c:function:: int CVodeSensReInit(void * cvode_mem, int ism, N_Vector * yS0)

   The routine :c:func:`CVodeSensReInit` reinitializes forward sensitivity computations.

   **Arguments:**
     * ``cvode_mem`` -- pointer to the CVODES memory block returned by :c:func:`CVodeCreate`.
     * ``ism`` --  forward sensitivity analysis!correction strategies a flag used to select the sensitivity solution method. Its value can be ``CV_SIMULTANEOUS`` , ``CV_STAGGERED`` , or ``CV_STAGGERED1``.
     * ``yS0`` -- a pointer to an array of ``Ns`` variables of type ``N_Vector`` containing the initial values of the sensitivities.

   **Return value:**
     * ``CV_SUCCESS`` -- The call to :c:func:`CVodeSensReInit` was successful.
     * ``CV_MEM_NULL`` -- The CVODES memory block was not initialized through a previous call to :c:func:`CVodeCreate`.
     * ``CV_NO_SENS`` -- Memory space for sensitivity integration was not allocated through a previous call to :c:func:`CVodeSensInit`.
     * ``CV_ILL_INPUT`` -- An input argument to :c:func:`CVodeSensReInit` has an illegal value.
     * ``CV_MEM_FAIL`` -- A memory allocation request has failed.

   **Notes:**
      All arguments of :c:func:`CVodeSensReInit` are the same as those of the
      functions  :c:func:`CVodeSensInit` and :c:func:`CVodeSensInit1`.  If an error
      occurred, :c:func:`CVodeSensReInit` also sends a message to the  error handler
      function.  :c:func:`CVodeSensReInit` potentially does some minimal memory
      allocation (for the  sensitivity absolute tolerance) and for arrays of
      counters used by the  ``CV_STAGGERED1`` method.

      .. warning::
         The value of the input argument ``ism`` must be compatible with  the
         type of the sensitivity ODE right-hand side function.  Thus  if the
         sensitivity module was initialized using :c:func:`CVodeSensInit`, then  it is
         illegal to pass ``ism`` = ``CV_STAGGERED1`` to :c:func:`CVodeSensReInit`.


To deallocate all forward sensitivity-related memory (allocated in a prior call
to :c:func:`CVodeSensInit` or :c:func:`CVodeSensInit1`), the user must call

.. c:function:: void CVodeSensFree(void * cvode_mem)

   The function :c:func:`CVodeSensFree` frees the memory allocated for forward
   sensitivity computations by a previous call to :c:func:`CVodeSensInit` or
   :c:func:`CVodeSensInit1`.

   **Arguments:**
     * ``cvode_mem`` -- pointer to the CVODES memory block returned by :c:func:`CVodeCreate`.

   **Return value:**
     * The function has no return value.

   **Notes:**
      In general, :c:func:`CVodeSensFree` need not be called by the user, as it is
      invoked automatically by :c:func:`CVodeFree`.

      After a call to :c:func:`CVodeSensFree`, forward sensitivity computations can be
      reactivated only by calling :c:func:`CVodeSensInit` or
      :c:func:`CVodeSensInit1` again.


To activate and deactivate forward sensitivity calculations for successive
CVODES runs, without having to allocate and deallocate memory, the following
function is provided:

.. c:function:: int CVodeSensToggleOff(void * cvode_mem)

   The function :c:func:`CVodeSensToggleOff` deactivates forward sensitivity
   calculations. It does not deallocate sensitivity-related memory.

   **Arguments:**
     * ``cvode_mem`` -- pointer to the memory previously returned by :c:func:`CVodeCreate`.

   **Return value:**
     * ``CV_SUCCESS`` -- :c:func:`CVodeSensToggleOff` was successful.
     * ``CV_MEM_NULL`` -- ``cvode_mem`` was ``NULL``.

   **Notes:**
      Since sensitivity-related memory is not deallocated, sensitivities can  be
      reactivated at a later time (using :c:func:`CVodeSensReInit`).


Forward sensitivity tolerance specification functions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

One of the following three functions must be called to specify the
integration tolerances for sensitivities. Note that this call must be made after
the call to :c:func:`CVodeSensInit` or :c:func:`CVodeSensInit1`.

.. c:function:: int CVodeSensSStolerances(void * cvode_mem, realtype reltolS, realtype* abstolS)

   The function :c:func:`CVodeSensSStolerances` specifies scalar relative and absolute
   tolerances.

   **Arguments:**
     * ``cvode_mem`` -- pointer to the CVODES memory block returned by :c:func:`CVodeCreate`.
     * ``reltolS`` -- is the scalar relative error tolerance.
     * ``abstolS`` -- is a pointer to an array of length ``Ns`` containing the scalar absolute error tolerances, one for each parameter.

   **Return value:**
     * ``CV_SUCCESS`` -- The call to ``CVodeSStolerances`` was successful.
     * ``CV_MEM_NULL`` -- The CVODES memory block was not initialized through a previous call to :c:func:`CVodeCreate`.
     * ``CV_NO_SENS`` -- The sensitivity allocation function :c:func:`CVodeSensInit` or :c:func:`CVodeSensInit1` has not been called.
     * ``CV_ILL_INPUT`` -- One of the input tolerances was negative.


.. c:function:: int CVodeSensSVtolerances(void * cvode_mem, realtype reltolS, N_Vector* abstolS)

   The function :c:func:`CVodeSensSVtolerances` specifies scalar relative tolerance
   and  vector absolute tolerances.

   **Arguments:**
     * ``cvode_mem`` -- pointer to the CVODES memory block returned by :c:func:`CVodeCreate`.
     * ``reltolS`` -- is the scalar relative error tolerance.
     * ``abstolS`` -- is an array of ``Ns`` variables of type ``N_Vector``. The ``N_Vector`` from ``abstolS[is]`` specifies the vector tolerances for ``is`` -th sensitivity.

   **Return value:**
     * ``CV_SUCCESS`` -- The call to ``CVodeSVtolerances`` was successful.
     * ``CV_MEM_NULL`` -- The CVODES memory block was not initialized through a previous call to :c:func:`CVodeCreate`.
     * ``CV_NO_SENS`` -- The allocation function for sensitivities has not been called.
     * ``CV_ILL_INPUT`` -- The relative error tolerance was negative or an absolute tolerance vector had a negative component.

   **Notes:**
      This choice of tolerances is important when the absolute error tolerance
      needs to  be different for each component of any vector ``yS[i]``.


.. c:function:: int CVodeSensEEtolerances(void * cvode_mem)

   When :c:func:`CVodeSensEEtolerances` is called, CVODES will estimate
   tolerances for  sensitivity variables based on the tolerances supplied for
   states variables  and the scaling factors :math:`\bar p`.

   **Arguments:**
     * ``cvode_mem`` -- pointer to the CVODES memory block returned by :c:func:`CVodeCreate`.

   **Return value:**
     * ``CV_SUCCESS`` -- The call to :c:func:`CVodeSensEEtolerances` was successful.
     * ``CV_MEM_NULL`` -- The CVODES memory block was not initialized through a previous call to :c:func:`CVodeCreate`.
     * ``CV_NO_SENS`` -- The sensitivity allocation function has not been called.


.. _CVODES.Usage.FSA.user_callable.nonlin_solv_init:

Forward sensitivity nonlinear solver interface functions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

As in the pure ODE case, when computing solution sensitivities using forward
sensitivitiy analysis CVODES uses the ``SUNNonlinearSolver`` implementation of
Newton’s method defined by the ``SUNNONLINSOL_NEWTON`` module (see
:numref:`SUNNonlinSol.Newton`) by default. To specify a different nonlinear
solver in CVODES, the user’s program must create a ``SUNNonlinearSolver`` object
by calling the appropriate constructor routine. The user must then attach the
``SUNNonlinearSolver`` object to CVODES by calling
:c:func:`CVodeSetNonlinearSolverSensSim` when using the ``CV_SIMULTANEOUS``
corrector option, or :c:func:`CVodeSetNonlinearSolver` and
:c:func:`CVodeSetNonlinearSolverSensStg` or
:c:func:`CVodeSetNonlinearSolverSensStg1` when using the ``CV_STAGGERED`` or
``CV_STAGGERED1`` corrector option respectively, as documented below.

When changing the nonlinear solver in CVODES, :c:func:`CVodeSetNonlinearSolver`
must be called after :c:func:`CVodeInit`; similarly
:c:func:`CVodeSetNonlinearSolverSensSim`, :c:func:`CVodeSetNonlinearSolverStg`,
and :c:func:`CVodeSetNonlinearSolverStg1` must be called after
:c:func:`CVodeSensInit`. If any calls to :c:func:`CVode` have been made, then CVODES
will need to be reinitialized by calling :c:func:`CVodeReInit` to ensure that
the nonlinear solver is initialized correctly before any subsequent calls to
:c:func:`CVode`.

The first argument passed to the routines
:c:func:`CVodeSetNonlinearSolverSensSim`,
:c:func:`CVodeSetNonlinearSolverSensStg`, and
:c:func:`CVodeSetNonlinearSolverSensStg1` is the CVODES memory pointer returned
by :c:func:`CVodeCreate` and the second argument is the ``SUNNonlinearSolver``
object to use for solving the nonlinear systems :eq:`CVODES_nonlinear` or
:eq:`CVODES_nonlinear_fixedpoint` A call to this function attaches the nonlinear solver
to the main CVODES integrator.


.. c:function:: int CVodeSetNonlinearSolverSensSim(void * cvode_mem, SUNNonlinearSolver NLS)

   The function :c:func:`CVodeSetNonLinearSolverSensSim` attaches a
   ``SUNNonlinearSolver``  object (``NLS``) to CVODES when using the
   ``CV_SIMULTANEOUS`` approach to  correct the state and sensitivity variables
   at the same time.

   **Arguments:**
     * ``cvode_mem`` -- pointer to the CVODES memory block.
     * ``NLS`` -- ``SUNNonlinearSolver`` object to use for solving nonlinear systems :eq:`CVODES_nonlinear` or :eq:`CVODES_nonlinear_fixedpoint`.

   **Return value:**
     * ``CV_SUCCESS`` -- The nonlinear solver was successfully attached.
     * ``CV_MEM_NULL`` -- The ``cvode_mem`` pointer is ``NULL``.
     * ``CV_ILL_INPUT`` -- The SUNNONLINSOL object is ``NULL``, does not implement the required nonlinear solver operations, is not of the correct type, or the residual function, convergence test function, or maximum number of nonlinear iterations could not be set.


.. c:function:: int CVodeSetNonlinearSolverSensStg(void * cvode_mem, SUNNonlinearSolver NLS)

   The function :c:func:`CVodeSetNonLinearSolverSensStg` attaches a
   ``SUNNonlinearSolver``  object (``NLS``) to CVODES when using the
   ``CV_STAGGERED`` approach to  correct all the sensitivity variables after the
   correction of the state  variables.

   **Arguments:**
     * ``cvode_mem`` -- pointer to the CVODES memory block.
     * ``NLS`` -- SUNNONLINSOL object to use for solving nonlinear systems.

   **Return value:**
     * ``CV_SUCCESS`` -- The nonlinear solver was successfully attached.
     * ``CV_MEM_NULL`` -- The ``cvode_mem`` pointer is ``NULL``.
     * ``CV_ILL_INPUT`` -- The SUNNONLINSOL object is ``NULL``, does not implement the required nonlinear solver operations, is not of the correct type, or the residual function, convergence test function, or maximum number of nonlinear iterations could not be set.

   **Notes:**
      This function only attaches the ``SUNNonlinearSolver`` object for
      correcting the  sensitivity variables. To attach a ``SUNNonlinearSolver``
      object for the state  variable correction use :c:func:`CVodeSetNonlinearSolver`.


.. c:function:: int CVodeSetNonlinearSolverSensStg1(void * cvode_mem, SUNNonlinearSolver NLS)

   The function :c:func:`CVodeSetNonLinearSolverSensStg1` attaches a
   ``SUNNonlinearSolver``  object (``NLS``) to CVODES when using the
   ``CV_STAGGERED1`` approach to  correct the sensitivity variables one at a
   time after the correction of the  state variables.

   **Arguments:**
     * ``cvode_mem`` -- pointer to the CVODES memory block.
     * ``NLS`` -- SUNNONLINSOL object to use for solving nonlinear systems.

   **Return value:**
     * ``CV_SUCCESS`` -- The nonlinear solver was successfully attached.
     * ``CV_MEM_NULL`` -- The ``cvode_mem`` pointer is ``NULL``.
     * ``CV_ILL_INPUT`` -- The SUNNONLINSOL object is ``NULL``, does not implement the required nonlinear solver operations, is not of the correct type, or the residual function, convergence test function, or maximum number of nonlinear iterations could not be set.

   **Notes:**
      This function only attaches the ``SUNNonlinearSolver`` object for
      correcting the  sensitivity variables. To attach a ``SUNNonlinearSolver``
      object for the state  variable correction use :c:func:`CVodeSetNonlinearSolver`.


CVODES solver function
^^^^^^^^^^^^^^^^^^^^^^

Even if forward sensitivity analysis was enabled, the call to the main solver function :c:func:`CVode` is exactly the same as
in :numref:`CVODES.Usage.SIM`. However, in this case the return value ``flag`` can also be one of the following:

- ``CV_SRHSFUNC_FAIL`` -- The sensitivity right-hand side function failed in an
  unrecoverable manner.

- ``CV_FIRST_SRHSFUNC_ERR`` -- The sensitivity right-hand side function failed
  at the first call.

- ``CV_REPTD_SRHSFUNC_ERR`` -- Convergence tests occurred too many times due to
  repeated recoverable errors in the sensitivity right-hand side
  function. This flag will also be returned if the sensitivity right-hand side
  function had repeated recoverable errors during the estimation of an initial
  step size.

- ``CV_UNREC_SRHSFUNC_ERR`` -- The sensitivity right-hand function had a
  recoverable error, but no recovery was possible. This failure mode is rare,
  as it can occur only if the sensitivity right-hand side function fails
  recoverably after an error test failed while at order one.


.. _CVODES.Usage.FSA.user_callable.sensi_get:

Forward sensitivity extraction functions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If forward sensitivity computations have been initialized by a call to
:c:func:`CVodeSensInit` or :c:func:`CVodeSensInit1`, or reinitialized by a call
to :c:func:`CVSensReInit`, then CVODES computes both a solution and
sensitivities at time ``t``. However, :c:func:`CVode` will still return only the
solution :math:`y` in ``yout``. Solution sensitivities can be obtained through
one of the following functions:


.. c:function:: int CVodeGetSens(void * cvode_mem, realtype * tret, N_Vector * yS)

   The function :c:func:`CVodeGetSens` returns the sensitivity solution vectors after
   a  successful return from :c:func:`CVode`.

   **Arguments:**
     * ``cvode_mem`` -- pointer to the memory previously allocated by :c:func:`CVodeInit`.
     * ``tret`` -- the time reached by the solver output.
     * ``yS`` -- array of computed forward sensitivity vectors. This vector array must be allocated by the user.

   **Return value:**
     * ``CV_SUCCESS`` -- :c:func:`CVodeGetSens` was successful.
     * ``CV_MEM_NULL`` -- ``cvode_mem`` was ``NULL``.
     * ``CV_NO_SENS`` -- Forward sensitivity analysis was not initialized.
     * ``CV_BAD_DKY`` -- ``yS`` is ``NULL``.

   **Notes:**
      Note that the argument ``tret`` is an output for this function. Its value
      will be  the same as that returned at the last :c:func:`CVode` call.


The function :c:func:`CVodeGetSensDky` computes the ``k``-th derivatives of the
interpolating polynomials for the sensitivity variables at time ``t``. This
function is called by :c:func:`CVodeGetSens` with ``k`` :math:`= 0`, but may also be
called directly by the user.


.. c:function:: int CVodeGetSensDky(void * cvode_mem, realtype t, int k, N_Vector * dkyS)

   The function :c:func:`CVodeGetSensDky` returns derivatives of the sensitivity
   solution  vectors after a successful return from :c:func:`CVode`.

   **Arguments:**
     * ``cvode_mem`` -- pointer to the memory previously allocated by :c:func:`CVodeInit`.
     * ``t`` -- specifies the time at which sensitivity information is requested. The time ``t`` must fall within the interval defined by the last successful step taken by CVODES.
     * ``k`` -- order of derivatives.
     * ``dkyS`` -- array of ``Ns`` vectors containing the derivatives on output. The space for ``dkyS`` must be allocated by the user.

   **Return value:**
     * ``CV_SUCCESS`` -- :c:func:`CVodeGetSensDky` succeeded.
     * ``CV_MEM_NULL`` -- ``cvode_mem`` was ``NULL``.
     * ``CV_NO_SENS`` -- Forward sensitivity analysis was not initialized.
     * ``CV_BAD_DKY`` -- One of the vectors ``dkyS`` is ``NULL``.
     * ``CV_BAD_K`` -- ``k`` is not in the range :math:`0, 1, ...,` ``qlast``.
     * ``CV_BAD_T`` -- The time ``t`` is not in the allowed range.


Forward sensitivity solution vectors can also be extracted separately for each
parameter in turn through the functions :c:func:`CVodeGetSens1` and
:c:func:`CVodeGetSensDky1`, defined as follows:


.. c:function:: int CVodeGetSens1(void * cvode_mem, realtype * tret, int is, N_Vector yS)

   The function :c:func:`CVodeGetSens1` returns the ``is``-th sensitivity solution
   vector  after a successful return from :c:func:`CVode`.

   **Arguments:**
     * ``cvode_mem`` -- pointer to the memory previously allocated by :c:func:`CVodeInit`.
     * ``tret`` -- the time reached by the solver output.
     * ``is`` -- specifies which sensitivity vector is to be returned :math:`0\le` ``is`` :math:`< N_s`.
     * ``yS`` -- the computed forward sensitivity vector. This vector array must be allocated by the user.

   **Return value:**
     * ``CV_SUCCESS`` -- :c:func:`CVodeGetSens1` was successful.
     * ``CV_MEM_NULL`` -- ``cvode_mem`` was ``NULL``.
     * ``CV_NO_SENS`` -- Forward sensitivity analysis was not initialized.
     * ``CV_BAD_IS`` -- The index ``is`` is not in the allowed range.
     * ``CV_BAD_DKY`` -- ``yS`` is ``NULL``.
     * ``CV_BAD_T`` -- The time ``t`` is not in the allowed range.

   **Notes:**
      Note that the argument ``tret`` is an output for this function. Its value
      will be  the same as that returned at the last :c:func:`CVode` call.


.. c:function:: int CVodeGetSensDky1(void * cvode_mem, realtype t, int k, int is, N_Vector dkyS)

   The function :c:func:`CVodeGetSensDky1` returns the ``k``-th derivative of the
   ``is``-th sensitivity solution vector after a successful return from
   :c:func:`CVode`.

   **Arguments:**
     * ``cvode_mem`` -- pointer to the memory previously allocated by :c:func:`CVodeInit`.
     * ``t`` -- specifies the time at which sensitivity information is requested. The time ``t`` must fall within the interval defined by the last successful step taken by CVODES.
     * ``k`` -- order of derivative.
     * ``is`` -- specifies the sensitivity derivative vector to be returned :math:`0\le` ``is`` :math:`< N_s`.
     * ``dkyS`` -- the vector containing the derivative. The space for ``dkyS`` must be allocated by the user.

   **Return value:**
     * ``CV_SUCCESS`` -- :c:func:`CVodeGetQuadDky1` succeeded.
     * ``CV_MEM_NULL`` -- The pointer to ``cvode_mem`` was ``NULL``.
     * ``CV_NO_SENS`` -- Forward sensitivity analysis was not initialized.
     * ``CV_BAD_DKY`` -- ``dkyS`` or one of the vectors ``dkyS[i]`` is ``NULL``.
     * ``CV_BAD_IS`` -- The index ``is`` is not in the allowed range.
     * ``CV_BAD_K`` -- ``k`` is not in the range :math:`0, 1, ...,` ``qlast``.
     * ``CV_BAD_T`` -- The time ``t`` is not in the allowed range.


.. _CVODES.Usage.FSA.user_callable.optional_inputs:

Optional inputs for forward sensitivity analysis
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Optional input variables that control the computation of sensitivities can be
changed from their default values through calls to ``CVodeSetSens*`` functions.
:numref:`CVODES.Usage.FSA.user_callable.optional_inputs.Table` lists all forward
sensitivity optional input functions in CVODES which are described in detail in
the remainder of this section.

We note that, on an error return, all of the optional input functions send an
error message to the error handler function. All error return values are
negative, so the test ``flag < 0`` will catch all errors. Finally, a call to a
``CVodeSetSens***`` function can be made from the user’s calling program at any
time and, if successful, takes effect immediately.

.. _CVODES.Usage.FSA.user_callable.optional_inputs.Table:
.. table:: Forward sensitivity optional inputs
   :align: center

   =================================== ==================================== ============
   **Optional input**                  **Routine name**                     **Default**
   =================================== ==================================== ============
   Sensitivity scaling factors         :c:func:`CVodeSetSensParams`         ``NULL``
   DQ approximation method             :c:func:`CVodeSetSensDQMethod`       centered/0.0
   Error control strategy              :c:func:`CVodeSetSensErrCon`         ``SUNFALSE``
   Maximum no. of nonlinear iterations :c:func:`CVodeSetSensMaxNonlinIters` 3
   =================================== ==================================== ============


.. c:function:: int CVodeSetSensParams(void * cvode_mem, realtype * p, realtype * pbar, int * plist)

   The function :c:func:`CVodeSetSensParams` specifies problem parameter information
   for sensitivity calculations.

   **Arguments:**
     * ``cvode_mem`` -- pointer to the CVODES memory block.
     * ``p`` -- a pointer to the array of real problem parameters used to
       evaluate :math:`f(t,y,p)`. If non-``NULL``, ``p`` must point to a field in
       the user's data structure ``user_data`` passed to the right-hand side
       function.
     * ``pbar`` -- an array of ``Ns`` positive scaling factors.
       If non-``NULL``, ``pbar`` must have all its components :math:`> 0.0`.
     * ``plist`` -- an array of ``Ns`` non-negative indices to specify
       which components ``p[i]`` to use in estimating the sensitivity equations.
       If non-``NULL``, ``plist`` must have all components :math:`\ge 0`.

   **Return value:**
     * ``CV_SUCCESS`` -- The optional value has been successfully set.
     * ``CV_MEM_NULL`` -- The ``cvode_mem`` pointer is ``NULL``.
     * ``CV_NO_SENS`` -- Forward sensitivity analysis was not initialized.
     * ``CV_ILL_INPUT`` -- An argument has an illegal value.

   **Notes:**
      .. warning::
         This function must be preceded by a call to :c:func:`CVodeSensInit` or
         :c:func:`CVodeSensInit1`.


.. c:function:: int CVodeSetSensDQMethod(void * cvode_mem, int DQtype, realtype DQrhomax)

   The function :c:func:`CVodeSetSensDQMethod` specifies the difference quotient
   strategy in  the case in which the right-hand side of the sensitivity
   equations are to  be computed by CVODES.

   **Arguments:**
     * ``cvode_mem`` -- pointer to the CVODES memory block.
     * ``DQtype`` -- specifies the difference quotient type. Its value can be ``CV_CENTERED`` or ``CV_FORWARD``.
     * ``DQrhomax`` -- positive value of the selection parameter used in deciding switching between a simultaneous or separate approximation of the two terms in the sensitivity right-hand side.

   **Return value:**
     * ``CV_SUCCESS`` -- The optional value has been successfully set.
     * ``CV_MEM_NULL`` -- The ``cvode_mem`` pointer is ``NULL``.
     * ``CV_ILL_INPUT`` -- An argument has an illegal value.

   **Notes:**
      If ``DQrhomax`` :math:`= 0.0`, then no switching is performed. The
      approximation is done  simultaneously using either centered or forward
      finite differences, depending on the  value of ``DQtype``.  For values of
      ``DQrhomax`` :math:`\ge 1.0`, the simultaneous  approximation is used
      whenever the estimated finite difference perturbations for  states and
      parameters are within a factor of ``DQrhomax``, and the separate
      approximation is used otherwise. Note that a value ``DQrhomax`` :math:`<1.0`
      will  effectively disable switching.   See :numref:`CVODES.Mathematics.FSA` for more details.
      The default value are ``DQtype == CV_CENTERED`` and
      ``DQrhomax=0.0``.


.. c:function:: int CVodeSetSensErrCon(void * cvode_mem, booleantype errconS)

   The function :c:func:`CVodeSetSensErrCon` specifies the error control  strategy for
   sensitivity variables.

   **Arguments:**
     * ``cvode_mem`` -- pointer to the CVODES memory block.
     * ``errconS`` -- specifies whether sensitivity variables are to be included ``SUNTRUE`` or not ``SUNFALSE`` in the error control mechanism.

   **Return value:**
     * ``CV_SUCCESS`` -- The optional value has been successfully set.
     * ``CV_MEM_NULL`` -- The ``cvode_mem`` pointer is ``NULL``.

   **Notes:**
      By default, ``errconS`` is set to ``SUNFALSE``.  If ``errconS = SUNTRUE``
      then both state variables and  sensitivity variables are included in the
      error tests.  If ``errconS = SUNFALSE`` then the sensitivity
      variables are excluded from the  error tests. Note that, in any event, all
      variables are considered in the convergence  tests.


.. c:function:: int CVodeSetSensMaxNonlinIters(void * cvode_mem, int maxcorS)

   The function :c:func:`CVodeSetSensMaxNonlinIters` specifies the maximum  number of
   nonlinear solver iterations for sensitivity variables per step.

   **Arguments:**
     * ``cvode_mem`` -- pointer to the CVODES memory block.
     * ``maxcorS`` -- maximum number of nonlinear solver iterations allowed per step :math:`> 0`.

   **Return value:**
     * ``CV_SUCCESS`` -- The optional value has been successfully set.
     * ``CV_MEM_NULL`` -- The ``cvode_mem`` pointer is ``NULL``.
     * ``CV_MEM_FAIL`` -- The SUNNONLINSOL module is ``NULL``.

   **Notes:**
      The default value is 3.


.. _CVODES.Usage.FSA.user_callable.optional_output:

Optional outputs for forward sensitivity analysis
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Optional output functions that return statistics and solver performance
information related to forward sensitivity computations are listed in
:numref:`CVODES.Usage.FSA.user_callable.optional_output.Table` and described in
detail in the remainder of this section.

.. _CVODES.Usage.FSA.user_callable.optional_output.Table:
.. table:: Forward sensitivity optional outputs
   :align: center

   ================================================== ================================================
   **Optional output**                                **Routine name**
   ================================================== ================================================
   No. of calls to sensitivity r.h.s. function        :c:func:`CVodeGetSensNumRhsEvals`
   No. of calls to r.h.s. function for sensitivity    :c:func:`CVodeGetNumRhsEvalsSens`
   No. of sensitivity local error test failures       :c:func:`CVodeGetSensNumErrTestFails`
   No. of calls to lin. solv. setup routine for sens. :c:func:`CVodeGetSensNumLinSolvSetups`
   Error weight vector for sensitivity variables      :c:func:`CVodeGetSensErrWeights`
   No. of sens. nonlinear solver iterations           :c:func:`CVodeGetSensNumNonlinSolvIters`
   No. of sens. convergence failures                  :c:func:`CVodeGetSensNumNonlinSolvConvFails`
   No. of staggered nonlinear solver iterations       :c:func:`CVodeGetStgrSensNumNonlinSolvIters`
   No. of staggered convergence failures              :c:func:`CVodeGetStgrSensNumNonlinSolvConvFails`
   ================================================== ================================================


.. c:function:: int CVodeGetSensNumRhsEvals(void * cvode_mem, long int nfSevals)

   The function :c:func:`CVodeGetSensNumRhsEvals` returns the number of calls to the
   sensitivity  right-hand side function.

   **Arguments:**
     * ``cvode_mem`` -- pointer to the CVODES memory block.
     * ``nfSevals`` -- number of calls to the sensitivity right-hand side function.

   **Return value:**
     * ``CV_SUCCESS`` -- The optional output value has been successfully set.
     * ``CV_MEM_NULL`` -- The ``cvode_mem`` pointer is ``NULL``.
     * ``CV_NO_SENS`` -- Forward sensitivity analysis was not initialized.

   **Notes:**
      In order to accommodate any of the three possible sensitivity solution
      methods,  the default internal  finite difference quotient functions
      evaluate the sensitivity right-hand sides  one at a time. Therefore,
      ``nfSevals`` will always be a multiple of the  number of sensitivity
      parameters (the same as the case in which the user supplies  a routine of
      type :c:type:`CVSensRhs1Fn`).


.. c:function:: int CVodeGetNumRhsEvalsSens(void * cvode_mem, long int nfevalsS)

   The function :c:func:`CVodeGetNumRhsEvalsSEns` returns the number of calls to the
   user's right-hand side function due to the internal finite difference
   approximation  of the sensitivity right-hand sides.

   **Arguments:**
     * ``cvode_mem`` -- pointer to the CVODES memory block.
     * ``nfevalsS`` -- number of calls to the user's ODE right-hand side function for the evaluation of sensitivity right-hand sides.

   **Return value:**
     * ``CV_SUCCESS`` -- The optional output value has been successfully set.
     * ``CV_MEM_NULL`` -- The ``cvode_mem`` pointer is ``NULL``.
     * ``CV_NO_SENS`` -- Forward sensitivity analysis was not initialized.

   **Notes:**
      This counter is incremented only if the internal finite difference
      approximation  routines are used for the evaluation of the sensitivity
      right-hand sides.


.. c:function:: int CVodeGetSensNumErrTestFails(void * cvode_mem, long int nSetfails)

   The function :c:func:`CVodeGetSensNumErrTestFails` returns the number of local
   error test failures for the sensitivity variables that have occurred.

   **Arguments:**
     * ``cvode_mem`` -- pointer to the CVODES memory block.
     * ``nSetfails`` -- number of error test failures.

   **Return value:**
     * ``CV_SUCCESS`` -- The optional output value has been successfully set.
     * ``CV_MEM_NULL`` -- The ``cvode_mem`` pointer is ``NULL``.
     * ``CV_NO_SENS`` -- Forward sensitivity analysis was not initialized.

   **Notes:**
      This counter is incremented only if the sensitivity variables have been
      included in the error test (see :c:func:`CVodeSetSensErrCon`).  Even in
      that case, this counter is not incremented if the
      ``ism = CV_SIMULTANEOUS``  sensitivity solution method has been used.


.. c:function:: int CVodeGetSensNumLinSolvSetups(void * cvode_mem, long int nlinsetupsS)

   The function :c:func:`CVodeGetSensNumLinSolvSetups` returns the number of calls  to the linear solver setup function due to forward sensitivity calculations.

   **Arguments:**
     * ``cvode_mem`` -- pointer to the CVODES memory block.
     * ``nlinsetupsS`` -- number of calls to the linear solver setup function.

   **Return value:**
     * ``CV_SUCCESS`` -- The optional output value has been successfully set.
     * ``CV_MEM_NULL`` -- The ``cvode_mem`` pointer is ``NULL``.
     * ``CV_NO_SENS`` -- Forward sensitivity analysis was not initialized.

   **Notes:**
      This counter is incremented only if a nonlinear solver requiring a linear
      solve has been used and if either the ``ism = CV_STAGGERED`` or the  ``ism
      = CV_STAGGERED1`` sensitivity solution method has been specified (see
      :numref:`CVODES.Usage.FSA.user_callable.sensi_malloc`).


.. c:function:: int CVodeGetSensStats(void* cvode_mem, long int* nfSevals, \
                long int* nfevalsS, long int* nSetfails, long int* nlinsetupsS)

   The function :c:func:`CVodeGetSensStats` returns all of the above
   sensitivity-related solver  statistics as a group.

   **Arguments:**
     * ``cvode_mem`` -- pointer to the CVODES memory block.
     * ``nfSevals`` -- number of calls to the sensitivity right-hand side function.
     * ``nfevalsS`` -- number of calls to the ODE right-hand side function for sensitivity evaluations.
     * ``nSetfails`` -- number of error test failures.
     * ``nlinsetupsS`` -- number of calls to the linear solver setup function.

   **Return value:**
     * ``CV_SUCCESS`` -- The optional output values have been successfully set.
     * ``CV_MEM_NULL`` -- The ``cvode_mem`` pointer is ``NULL``.
     * ``CV_NO_SENS`` -- Forward sensitivity analysis was not initialized.


.. c:function:: int CVodeGetSensErrWeights(void * cvode_mem, N_Vector * eSweight)

   The function :c:func:`CVodeGetSensErrWeights` returns the sensitivity error weight
   vectors at the current time. These are the reciprocals of the :math:`W_i` of
   :eq:`CVODES_errwt` for the sensitivity variables.

   **Arguments:**
     * ``cvode_mem`` -- pointer to the CVODES memory block.
     * ``eSweight`` -- pointer to the array of error weight vectors.

   **Return value:**
     * ``CV_SUCCESS`` -- The optional output value has been successfully set.
     * ``CV_MEM_NULL`` -- The ``cvode_mem`` pointer is ``NULL``.
     * ``CV_NO_SENS`` -- Forward sensitivity analysis was not initialized.

   **Notes:**
      The user must allocate memory for ``eweightS``.


.. c:function:: int CVodeGetSensNumNonlinSolvIters(void * cvode_mem, long int nSniters)

   The function :c:func:`CVodeGetSensNumNonlinSolvIters` returns the  number of
   nonlinear iterations performed for  sensitivity calculations.

   **Arguments:**
     * ``cvode_mem`` -- pointer to the CVODES memory block.
     * ``nSniters`` -- number of nonlinear iterations performed.

   **Return value:**
     * ``CV_SUCCESS`` -- The optional output value has been successfully set.
     * ``CV_MEM_NULL`` -- The ``cvode_mem`` pointer is ``NULL``.
     * ``CV_NO_SENS`` -- Forward sensitivity analysis was not initialized.
     * ``CV_MEM_FAIL`` -- The SUNNONLINSOL module is ``NULL``.

   **Notes:**
      This counter is incremented only if ``ism`` was ``CV_STAGGERED`` or
      ``CV_STAGGERED1`` (see
      :numref:`CVODES.Usage.FSA.user_callable.sensi_malloc`).  In the
      ``CV_STAGGERED1`` case, the value of ``nSniters`` is the sum of  the
      number of nonlinear iterations performed for each sensitivity equation.
      These individual counters can be obtained through a call to
      :c:func:`CVodeGetStgrSensNumNonlinSolvIters` (see below).


.. c:function:: int CVodeGetSensNumNonlinSolvConvFails(void * cvode_mem, long int nSncfails)

   The function :c:func:`CVodeGetSensNumNonlinSolvConvFails` returns the  number of
   nonlinear convergence failures that have occurred for  sensitivity
   calculations.

   **Arguments:**
     * ``cvode_mem`` -- pointer to the CVODES memory block.
     * ``nSncfails`` -- number of nonlinear convergence failures.

   **Return value:**
     * ``CV_SUCCESS`` -- The optional output value has been successfully set.
     * ``CV_MEM_NULL`` -- The ``cvode_mem`` pointer is ``NULL``.
     * ``CV_NO_SENS`` -- Forward sensitivity analysis was not initialized.

   **Notes:**
      This counter is incremented only if ``ism`` was ``CV_STAGGERED`` or
      ``CV_STAGGERED1``.  In the ``CV_STAGGERED1`` case, the value of
      ``nSncfails`` is the sum of  the number of nonlinear convergence failures
      that occurred for each sensitivity equation.  These individual counters
      can be obtained through a call to  :c:func:`CVodeGetStgrSensNumNonlinConvFails`
      (see below).


.. c:function:: int CVodeGetSensNonlinSolvStats(void * cvode_mem, long int nSniters, long int nSncfails)

   The function :c:func:`CVodeGetSensNonlinSolvStats` returns the sensitivity-related
   nonlinear solver statistics as a group.

   **Arguments:**
     * ``cvode_mem`` -- pointer to the CVODES memory block.
     * ``nSniters`` -- number of nonlinear iterations performed.
     * ``nSncfails`` -- number of nonlinear convergence failures.

   **Return value:**
     * ``CV_SUCCESS`` -- The optional output values have been successfully set.
     * ``CV_MEM_NULL`` -- The ``cvode_mem`` pointer is ``NULL``.
     * ``CV_NO_SENS`` -- Forward sensitivity analysis was not initialized.
     * ``CV_MEM_FAIL`` -- The SUNNONLINSOL module is ``NULL``.


.. c:function:: int CVodeGetStgrSensNumNonlinSolvIters(void * cvode_mem, long int * nSTGR1niters)

   The function :c:func:`CVodeGetStgrSensNumNonlinSolvIters` returns the  number of
   nonlinear iterations performed for  each sensitivity equation separately, in
   the ``CV_STAGGERED1`` case.

   **Arguments:**
     * ``cvode_mem`` -- pointer to the CVODES memory block.
     * ``nSTGR1niters`` -- an array of dimension ``Ns`` which will be set with the number of nonlinear iterations performed for each sensitivity system individually.

   **Return value:**
     * ``CV_SUCCESS`` -- The optional output value has been successfully set.
     * ``CV_MEM_NULL`` -- The ``cvode_mem`` pointer is ``NULL``.
     * ``CV_NO_SENS`` -- Forward sensitivity analysis was not initialized.

   **Notes:**
      .. warning::
         The user must allocate space for ``nSTGR1niters``.


.. c:function:: int CVodeGetStgrSensNumNonlinSolvConvFails(void * cvode_mem, long int * nSTGR1ncfails)

   The function :c:func:`CVodeGetStgrSensNumNonlinSolvConvFails` returns the  number
   of nonlinear convergence failures that have occurred for  each sensitivity
   equation separately, in the ``CV_STAGGERED1`` case.

   **Arguments:**
     * ``cvode_mem`` -- pointer to the CVODES memory block.
     * ``nSTGR1ncfails`` -- an array of dimension ``Ns`` which will be set with the number of nonlinear convergence failures for each sensitivity system individually.

   **Return value:**
     * ``CV_SUCCESS`` -- The optional output value has been successfully set.
     * ``CV_MEM_NULL`` -- The ``cvode_mem`` pointer is ``NULL``.
     * ``CV_NO_SENS`` -- Forward sensitivity analysis was not initialized.

   **Notes:**
      .. warning::
         The user must allocate space for ``nSTGR1ncfails``.


.. c:function:: int CVodeGetStgrSensNonlinSolvStats(void * cvode_mem, long int * nSTRG1niterslong, int * nSTGR1ncfails)

   The function :c:func:`CVodeGetStgrSensNonlinSolvStats` returns the  number of
   nonlinear iterations and convergence failures that have occurred for  each
   sensitivity equation separately, in the ``CV_STAGGERED1`` case.

   **Arguments:**
     * ``cvode_mem`` -- pointer to the CVODES memory block.
     * ``nSTGR1niters`` -- an array of dimension ``Ns`` which will be set with the number of nonlinear iterations performed for each sensitivity system individually.
     * ``nSTGR1ncfails`` -- an array of dimension ``Ns`` which will be set with the number of nonlinear convergence failures for each sensitivity system individually.

   **Return value:**
     * ``CV_SUCCESS`` -- The optional output values have been successfully set.
     * ``CV_MEM_NULL`` -- The ``cvode_mem`` pointer is ``NULL``.
     * ``CV_NO_SENS`` -- Forward sensitivity analysis was not initialized.
     * ``CV_MEM_FAIL`` -- The SUNNONLINSOL module is ``NULL``.


.. _CVODES.Usage.FSA.user_supplied:

User-supplied routines for forward sensitivity analysis
-------------------------------------------------------

In addition to the required and optional user-supplied routines described in
:numref:`CVODES.Usage.SIM.user_supplied`, when using CVODES for forward
sensitivity analysis, the user has the option of providing a routine that
calculates the right-hand side of the sensitivity equations :eq:`CVODES_sens_eqns`.

By default, CVODES uses difference quotient approximation routines for the
right-hand sides of the sensitivity equations. However, CVODES allows the option
for user-defined sensitivity right-hand side routines (which also provides a
mechanism for interfacing CVODES to routines generated by automatic
differentiation).

Sensitivity equations right-hand side (all at once)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If the ``CV_SIMULTANEOUS`` or ``CV_STAGGERED`` approach was selected in the call
to :c:func:`CVodeSensInit` or :c:func:`CVodeSensInit1`, the user may provide the
right-hand sides of the sensitivity equations :eq:`CVODES_sens_eqns`, for all
sensitivity parameters at once, through a function of type :c:type:`CVSensRhsFn`
defined by:

.. c:type:: int (*CVSensRhsFn)(int Ns, realtype t, N_Vector y, N_Vector ydot, N_Vector *yS, N_Vector *ySdot, void *user_data, N_Vector tmp1, N_Vector tmp2)

   This function computes the sensitivity right-hand side for all sensitivity
   equations at once.  It must compute the vectors
   :math:`\dfrac{\partial f}{\partial y} s_i(t) + \dfrac{\partial f}{\partial p_i}`
   and store them in ``ySdot[i]``.

   **Arguments:**
     * ``Ns`` -- is the number of sensitivities.
     * ``t`` -- is the current value of the independent variable.
     * ``y`` -- is the current value of the state vector, :math:`y(t)` .
     * ``ydot`` -- is the current value of the right-hand side of the state equations.
     * ``yS`` -- contains the current values of the sensitivity vectors.
     * ``ySdot`` -- is the output of :c:type:`CVSensRhsFn` . On exit it must contain    the sensitivity right-hand side vectors.
     * ``user_data`` -- is a pointer to user data, the same as the ``user_data`` parameter passed to :c:func:`CVodeSetUserData` .
     * ``tmp1``, ``tmp2`` -- are ``N_Vectors`` of length :math:`N` which can be used as temporary storage.

   **Return value:**
      A :c:type:`CVSensRhsFn` should return 0 if successful, a positive value if a recoverable
      error occurred (in which case CVODES will attempt to correct), or a negative
      value if it failed unrecoverably (in which case the integration is halted and
      ``CV_SRHSFUNC_FAIL`` is returned).

   **Notes:**
       Allocation of memory for ``ySdot`` is handled within CVODES.  There are
       two situations in which recovery is not possible even if
       :c:type:`CVSensRhsFn`  function returns a recoverable error flag.  One is
       when this  occurs at the very first call to the :c:type:`CVSensRhsFn` (in
       which case CVODES returns  ``CV_FIRST_SRHSFUNC_ERR``).  The other is when
       a recoverable error is reported  by :c:type:`CVSensRhsFn` after an error
       test failure, while the linear multistep method  order is equal to 1 (in
       which case CVODES returns ``CV_UNREC_SRHSFUNC_ERR``).

       .. warning::
          A sensitivity right-hand side function of type :c:type:`CVSensRhsFn`
          is not compatible with the ``CV_STAGGERED1`` approach.


Sensitivity equations right-hand side (one at a time)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Alternatively, the user may provide the sensitivity right-hand sides, one
sensitivity parameter at a time, through a function of type
:c:type:`CVSensRhs1Fn`. Note that a sensitivity right-hand side function of type
:c:type:`CVSensRhs1Fn` is compatible with any valid value of the argument
``ism`` to :c:func:`CVodeSensInit` and :c:func:`CVodeSensInit1`, and is
*required* if ``ism = CV_STAGGERED1`` in the call to
:c:func:`CVodeSensInit1`. The type :c:type:`CVSensRhs1Fn` is defined by


.. c:type:: int (*CVSensRhs1Fn)(int Ns, realtype t, N_Vector y, N_Vector ydot, int iS, N_Vector yS, N_Vector ySdot, void *user_data, N_Vector tmp1, N_Vector tmp2)

   This function computes the sensitivity right-hand side for one sensitivity
   equation at a time.  It must compute the vector
   :math:`(\frac{\partial f}{\partial y}) s_i(t) + (\frac{\partial f}{\partial p_i})`
   for :math:`i` = ``iS`` and  store it in ``ySdot``.

   **Arguments:**
     * ``Ns`` -- is the number of sensitivities.
     * ``t`` -- is the current value of the independent variable.
     * ``y`` -- is the current value of the state vector, :math:`y(t)` .
     * ``ydot`` -- is the current value of the right-hand side of the state equations.
     * ``iS`` -- is the index of the parameter for which the sensitivity right-hand    side must be computed :math:`(0 \leq`  ``iS``  :math:`<`  ``Ns``).
     * ``yS`` -- contains the current value of the ``iS`` -th sensitivity vector.
     * ``ySdot`` -- is the output of :c:type:`CVSensRhs1Fn` . On exit it must contain    the ``iS`` -th sensitivity right-hand side vector.
     * ``user_data`` -- is a pointer to user data, the same as the ``user_data`` parameter passed to :c:func:`CVodeSetUserData` .
     * ``tmp1``, ``tmp2`` -- are ``N_Vectors`` of length :math:`N` which can be used as temporary storage.

   **Return value:**
      A :c:type:`CVSensRhs1Fn` should return 0 if successful, a positive value if a recoverable
      error occurred (in which case CVODES will attempt to correct), or a negative
      value if it failed unrecoverably (in which case the integration is halted and
      ``CV_SRHSFUNC_FAIL`` is returned).

   **Notes:**
      Allocation of memory for ``ySdot`` is handled within CVODES.  There
      are two situations in which recovery is not possible even if
      :c:type:`CVSensRhs1Fn`  function returns a recoverable error flag.  One is when
      this occurs  at the very first call to the :c:type:`CVSensRhs1Fn` (in which case
      CVODES returns  ``CV_FIRST_SRHSFUNC_ERR``).  The other is when a
      recoverable error is reported  by :c:type:`CVSensRhs1Fn` after an error test
      failure, while the linear multistep method  order  equal to 1 (in which
      case CVODES returns ``CV_UNREC_SRHSFUNC_ERR``).


.. _CVODES.Usage.FSA.quad:

Integration of quadrature equations depending on forward sensitivities
----------------------------------------------------------------------

CVODES provides support for integration of quadrature equations that depends not
only on the state variables but also on forward sensitivities.

The following is an overview of the sequence of calls in a user’s main program
in this situation. Steps that are unchanged from the skeleton program presented
in :numref:`CVODES.Usage.FSA.skeleton_sim` are left unbolded.

#. Initialize parallel or multi-threaded environment, if appropriate

#. Create ``SUNContext`` object by calling :c:func:`SUNContext_Create`

#. Set problem dimensions etc.

#. Set vectors of initial values

#. Create CVODES object

#. Initialize CVODES solver

#. Specify integration tolerances

#. Create matrix object

#. Create linear solver object

#. Set linear solver optional inputs

#. Attach linear solver module

#. Set optional inputs

#. Create nonlinear solver object

#. Attach nonlinear solver module

#. Set nonlinear solver optional inputs

#. Initialize sensitivity-independent quadrature problem

#. Define the sensitivity problem

#. Set sensitivity initial conditions

#. Activate sensitivity calculations

#. Set sensitivity tolerances

#. Set sensitivity analysis optional inputs

#. Create sensitivity nonlinear solver object

#. Attach the sensitvity nonlinear solver module

#. Set sensitivity nonlinear solver optional inputs

#. **Set vector of initial values for quadrature variables**

#. Typically, the quadrature variables should be initialized to :math:`0`.

#. **Initialize sensitivity-dependent quadrature integration**

#. Call :c:func:`CVodeQuadSensInit` to specify the quadrature equation
   right-hand side function and to allocate internal memory related to
   quadrature integration.

#. **Set optional inputs for sensitivity-dependent quadrature integration**

#. Call :c:func:`CVodeSetQuadSensErrCon` to indicate whether or not quadrature variables should be used in the step size
   control mechanism. If so, one of the ``CVodeQuadSens*tolerances`` functions must be called to specify the integration
   tolerances for quadrature variables.

#. Advance solution in time

#. **Extract sensitivity-dependent quadrature variables**

#. Call :c:func:`CVodeGetQuadSens`, :c:func:`CVodeGetQuadSens1`, :c:func:`CVodeGetQuadSensDky` or :c:func:`CVodeGetQuadSensDky1` to obtain the
   values of the quadrature variables or their derivatives at the current time.

#. Get optional outputs

#. Extract sensitivity solution

#. **Get sensitivity-dependent quadrature optional outputs**

#. Call ``CVodeGetQuadSens*`` functions to obtain desired optional output related to the integration of
   sensitivity-dependent quadratures.

#. Deallocate memory for solutions vector

#. Deallocate memory for sensitivity vectors

#. **Deallocate memory for sensitivity-dependent quadrature variables**

#. **Free solver memory**

#. Free nonlinear solver memory

#. Free vector specification memory

#. Free linear solver and matrix memory

#. Free ``SUNContext`` object with a call to :c:func:`SUNContext_Free`

#. Finalize MPI, if used


.. _CVODES.Usage.FSA.quad.quad_sens_init:

Sensitivity-dependent quadrature initialization and deallocation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The function :c:func:`CVodeQuadSensInit` activates integration of quadrature equations
depending on sensitivities and allocates internal memory related to these
calculations. If ``rhsQS`` is input as ``NULL``, then CVODES uses an internal
function that computes difference quotient approximations to the functions
:math:`\bar{q}_i = q_y s_i + q_{p_i}`, in the notation of :eq:`CVODES_QUAD`. The form
of the call to this function is as follows:


.. c:function:: int CVodeQuadSensInit(void * cvode_mem, CVQuadSensRhsFn rhsQS, N_Vector * yQS0)

   The function :c:func:`CVodeQuadSensInit` provides required problem
   specifications,  allocates internal memory, and initializes quadrature
   integration.

   **Arguments:**
     * ``cvode_mem`` -- pointer to the CVODES memory block returned by :c:func:`CVodeCreate`.
     * ``rhsQS`` -- is the function which computes :math:`f_{QS}` , the right-hand side of the sensitivity-dependent quadrature..
     * ``yQS0`` -- contains the initial values of sensitivity-dependent quadratures.

   **Return value:**
     * ``CV_SUCCESS`` -- The call to :c:func:`CVodeQuadSensInit` was successful.
     * ``CVODE_MEM_NULL`` -- The CVODES memory was not initialized by a prior call to :c:func:`CVodeCreate`.
     * ``CVODE_MEM_FAIL`` -- A memory allocation request failed.
     * ``CV_NO_SENS`` -- The sensitivities were not initialized by a prior call to :c:func:`CVodeSensInit` or :c:func:`CVodeSensInit1`.
     * ``CV_ILL_INPUT`` -- The parameter ``yQS0`` is ``NULL``.

   **Notes:**
      .. warning::
          Before calling :c:func:`CVodeQuadSensInit`, the user must enable the
          sensitivites  by calling :c:func:`CVodeSensInit` or :c:func:`CVodeSensInit1`.  If
          an error occurred, :c:func:`CVodeQuadSensInit` also sends an error
          message to the  error handler function.


.. c:function:: int CVodeQuadSensReInit(void * cvode_mem, N_Vector * yQS0)

   The function :c:func:`CVodeQuadSensReInit` provides required problem
   specifications  and reinitializes the sensitivity-dependent quadrature
   integration.

   **Arguments:**
     * ``cvode_mem`` -- pointer to the CVODES memory block.
     * ``yQS0`` -- contains the initial values of sensitivity-dependent quadratures.

   **Return value:**
     * ``CV_SUCCESS`` -- The call to :c:func:`CVodeQuadSensReInit` was successful.
     * ``CVODE_MEM_NULL`` -- The CVODES memory was not initialized by a prior call to :c:func:`CVodeCreate`.
     * ``CV_NO_SENS`` -- Memory space for the sensitivity calculation was not allocated by a prior call to :c:func:`CVodeSensInit` or :c:func:`CVodeSensInit1`.
     * ``CV_NO_QUADSENS`` -- Memory space for the sensitivity quadratures integration was not allocated by a prior call to :c:func:`CVodeQuadSensInit`.
     * ``CV_ILL_INPUT`` -- The parameter ``yQS0`` is ``NULL``.

   **Notes:**
      If an error occurred, :c:func:`CVodeQuadSensReInit` also sends an error
      message to the  error handler function.


.. c:function:: void CVodeQuadSensFree(void* cvode_mem)

   The function :c:func:`CVodeQuadSensFree` frees the memory allocated for
   sensitivity quadrature integration.

   **Arguments:**
     * ``cvode_mem`` -- pointer to the CVODE memory block.

   **Return value:**
      There is no return value.

   **Notes:**
      In general, :c:func:`CVodeQuadSensFree` need not be called by the user as it
      is called automatically by :c:func:`CVodeFree`.


CVODES solver function
^^^^^^^^^^^^^^^^^^^^^^

Even if quadrature integration was enabled, the call to the main solver function
:c:func:`CVode` is exactly the same as in :numref:`CVODES.Usage.SIM`. However,
in this case the return value ``flag`` can also be one of the following:

- ``CV_QSRHSFUNC_ERR`` -- The sensitivity quadrature right-hand side
   function failed in an unrecoverable manner.

- ``CV_FIRST_QSRHSFUNC_ERR`` -- The sensitivity quadrature right-hand side
   function failed at the first call.

- ``CV_REPTD_QSRHSFUNC_ERR`` -- Convergence test failures occurred too many
  times due to repeated recoverable errors in the quadrature right-hand side
  function. This flag will also be returned if the quadrature right-hand side
  function had repeated recoverable errors during the estimation of an initial
  step size (assuming the sensitivity quadrature variables are included in the
  error tests).


.. _CVODES.Usage.FSA.quad.quad_sens_get:

Sensitivity-dependent quadrature extraction functions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If sensitivity-dependent quadratures have been initialized by a call to
:c:func:`CVodeQuadSensInit`, or reinitialized by a call to
:c:func:`CVodeQuadSensReInit`, then CVODES computes a solution, sensitivity
vectors, and quadratures depending on sensitivities at time ``t``. However,
:c:func:`CVode` will still return only the solution :math:`y`.
Sensitivity-dependent quadratures can be obtained using one of the following
functions:

.. c:function:: int CVodeGetQuadSens(void * cvode_mem, realtype tret, N_Vector * yQS)

   The function :c:func:`CVodeGetQuadSens` returns the quadrature sensitivities
   solution vectors after a successful return from :c:func:`CVode`.

   **Arguments:**
     * ``cvode_mem`` -- pointer to the memory previously allocated by :c:func:`CVodeInit`.
     * ``tret`` -- the time reached by the solver output.
     * ``yQS`` -- array of ``Ns`` computed sensitivity-dependent quadrature vectors. This vector array must be allocated by the user.

   **Return value:**
     * ``CV_SUCCESS`` -- :c:func:`CVodeGetQuadSens` was successful.
     * ``CVODE_MEM_NULL`` -- ``cvode_mem`` was ``NULL``.
     * ``CV_NO_SENS`` -- Sensitivities were not activated.
     * ``CV_NO_QUADSENS`` -- Quadratures depending on the sensitivities were not activated.
     * ``CV_BAD_DKY`` -- ``yQS`` or one of the ``yQS[i]`` is ``NULL``.


The function :c:func:`CVodeGetQuadSensDky` computes the ``k``-th derivatives of
the interpolating polynomials for the sensitivity-dependent quadrature variables
at time ``t``. This function is called by :c:func:`CVodeGetQuadSens` with ``k =
0``, but may also be called directly by the user.

.. c:function:: int CVodeGetQuadSensDky(void* cvode_mem, realtype t, int k, N_Vector* dkyQS)

   The function :c:func:`CVodeGetQuadSensDky` returns derivatives of the
   quadrature sensitivities  solution vectors after a successful return from
   :c:func:`CVode`.

   **Arguments:**
     * ``cvode_mem`` -- pointer to the memory previously allocated by :c:func:`CVodeInit`.
     * ``t`` -- the time at which information is requested. The time ``t`` must fall within the interval defined by the last successful step taken by CVODES.
     * ``k`` -- order of the requested derivative.
     * ``dkyQS`` -- array of ``Ns`` the vector containing the derivatives on output. This vector array must be allocated by the user.

   **Return value:**
     * ``CV_SUCCESS`` -- :c:func:`CVodeGetQuadSensDky` succeeded.
     * ``CVODE_MEM_NULL`` -- The pointer to ``cvode_mem`` was ``NULL``.
     * ``CV_NO_SENS`` -- Sensitivities were not activated.
     * ``CV_NO_QUADSENS`` -- Quadratures depending on the sensitivities were not activated.
     * ``CV_BAD_DKY`` -- ``dkyQS`` or one of the vectors ``dkyQS[i]`` is ``NULL``.
     * ``CV_BAD_K`` -- ``k`` is not in the range :math:`0, 1, ...,` ``qlast``.
     * ``CV_BAD_T`` -- The time ``t`` is not in the allowed range.


Quadrature sensitivity solution vectors can also be extracted separately for
each parameter in turn through the functions :c:func:`CVodeGetQuadSens1` and
:c:func:`CVodeGetQuadSensDky1`, defined as follows:

.. c:function:: int CVodeGetQuadSens1(void * cvode_mem, realtype tret, int is, N_Vector yQS)

   The function :c:func:`CVodeGetQuadSens1` returns the ``is``-th sensitivity
   of quadratures after a successful return from :c:func:`CVode`.

   **Arguments:**
     * ``cvode_mem`` -- pointer to the memory previously allocated by :c:func:`CVodeInit`.
     * ``tret`` -- the time reached by the solver output.
     * ``is`` -- specifies which sensitivity vector is to be returned :math:`0 \le` ``is`` :math:`< N_s`.
     * ``yQS`` -- the computed sensitivity-dependent quadrature vector. This vector array must be allocated by the user.

   **Return value:**
     * ``CV_SUCCESS`` -- :c:func:`CVodeGetQuadSens1` was successful.
     * ``CVODE_MEM_NULL`` -- ``cvode_mem`` was ``NULL``.
     * ``CV_NO_SENS`` -- Forward sensitivity analysis was not initialized.
     * ``CV_NO_QUADSENS`` -- Quadratures depending on the sensitivities were not activated.
     * ``CV_BAD_IS`` -- The index ``is`` is not in the allowed range.
     * ``CV_BAD_DKY`` -- ``yQS`` is ``NULL``.


.. c:function:: int CVodeGetQuadSensDky1(void * cvode_mem, realtype t, int k, int is, N_Vector dkyQS)

   The function :c:func:`CVodeGetQuadSensDky1` returns the ``k``-th derivative
   of the  ``is``-th sensitivity solution vector after a successful  return from
   :c:func:`CVode`.

   **Arguments:**
     * ``cvode_mem`` -- pointer to the memory previously allocated by :c:func:`CVodeInit`.
     * ``t`` -- specifies the time at which sensitivity information is requested. The time ``t`` must fall within the interval defined by the last successful step taken by CVODES.
     * ``k`` -- order of derivative.
     * ``is`` -- specifies the sensitivity derivative vector to be returned :math:`0\le` ``is`` :math:`< N_s`.
     * ``dkyQS`` -- the vector containing the derivative on output. The space for ``dkyQS`` must be allocated by the user.

   **Return value:**
     * ``CV_SUCCESS`` -- :c:func:`CVodeGetQuadDky1` succeeded.
     * ``CVODE_MEM_NULL`` -- ``cvode_mem`` was ``NULL``.
     * ``CV_NO_SENS`` -- Forward sensitivity analysis was not initialized.
     * ``CV_NO_QUADSENS`` -- Quadratures depending on the sensitivities were not activated.
     * ``CV_BAD_DKY`` -- ``dkyQS`` is ``NULL``.
     * ``CV_BAD_IS`` -- The index ``is`` is not in the allowed range.
     * ``CV_BAD_K`` -- ``k`` is not in the range :math:`0, 1, ...,` ``qlast``.
     * ``CV_BAD_T`` -- The time ``t`` is not in the allowed range.


.. _CVODES.Usage.FSA.quad.quad_sens_optional_input:

Optional inputs for sensitivity-dependent quadrature integration
----------------------------------------------------------------

CVODES provides the following optional input functions to control the
integration of sensitivity-dependent quadrature equations.


.. c:function:: int CVodeSetQuadSensErrCon(void * cvode_mem, booleantype errconQS)

   The function :c:func:`CVodeSetQuadSensErrCon` specifies whether or not the
   quadrature variables are to be used in the step size control  mechanism. If
   they are, the user must call one of the functions
   :c:func:`CVodeQuadSensSStolerances`, :c:func:`CVodeQuadSensSVtolerances`, or
   :c:func:`CVodeQuadSensEEtolerances` to specify the integration tolerances for
   the quadrature variables.

   **Arguments:**
     * ``cvode_mem`` -- pointer to the CVODES memory block.
     * ``errconQS`` -- specifies whether sensitivity quadrature variables are to be included ``SUNTRUE`` or not ``SUNFALSE`` in the error control mechanism.

   **Return value:**
     * ``CV_SUCCESS`` -- The optional value has been successfully set.
     * ``CVODE_MEM_NULL`` -- ``cvode_mem`` is ``NULL``.
     * ``CV_NO_SENS`` -- Sensitivities were not activated.
     * ``CV_NO_QUADSENS`` -- Quadratures depending on the sensitivities were not activated.

   **Notes:**
      By default, ``errconQS`` is set to ``SUNFALSE``.

      .. warning::
         It is illegal to call :c:func:`CVodeSetQuadSensErrCon` before a call
         to :c:func:`CVodeQuadSensInit`.


.. c:function:: int CVodeQuadSensSStolerances(void * cvode_mem, realtype reltolQS, realtype* abstolQS)

   The function :c:func:`CVodeQuadSensSStolerances` specifies scalar relative
   and absolute  tolerances.

   **Arguments:**
     * ``cvode_mem`` -- pointer to the CVODES memory block.
     * ``reltolQS`` --  tolerances is the scalar relative error tolerance.
     * ``abstolQS`` -- is a pointer to an array containing the ``Ns`` scalar absolute error tolerances.

   **Return value:**
     * ``CV_SUCCESS`` -- The optional value has been successfully set.
     * ``CVODE_MEM_NULL`` -- The ``cvode_mem`` pointer is ``NULL``.
     * ``CV_NO_SENS`` -- Sensitivities were not activated.
     * ``CV_NO_QUADSENS`` -- Quadratures depending on the sensitivities were not activated.
     * ``CV_ILL_INPUT`` -- One of the input tolerances was negative.

.. c:function:: int CVodeQuadSensSVtolerances(void * cvode_mem, realtype reltolQS, N_Vector* abstolQS)

   The function :c:func:`CVodeQuadSensSVtolerances` specifies scalar relative
   and  vector absolute tolerances.

   **Arguments:**
     * ``cvode_mem`` -- pointer to the CVODES memory block.
     * ``reltolQS`` --  tolerances is the scalar relative error tolerance.
     * ``abstolQS`` -- is an array of ``Ns`` variables of type ``N_Vector``. The ``N_Vector`` ``abstolS[is]`` specifies the vector tolerances for ``is`` -th quadrature sensitivity.

   **Return value:**
     * ``CV_SUCCESS`` -- The optional value has been successfully set.
     * ``CV_NO_QUAD`` -- Quadrature integration was not initialized.
     * ``CVODE_MEM_NULL`` -- The ``cvode_mem`` pointer is ``NULL``.
     * ``CV_NO_SENS`` -- Sensitivities were not activated.
     * ``CV_NO_QUADSENS`` -- Quadratures depending on the sensitivities were not activated.
     * ``CV_ILL_INPUT`` -- One of the input tolerances was negative.


.. c:function:: int CVodeQuadSensEEtolerances(void * cvode_mem)

   A call to the function :c:func:`CVodeQuadSensEEtolerances` specifies that the
   tolerances for the sensitivity-dependent quadratures should be estimated from
   those provided for the pure quadrature variables.

   **Arguments:**
     * ``cvode_mem`` -- pointer to the CVODES memory block.

   **Return value:**
     * ``CV_SUCCESS`` -- The optional value has been successfully set.
     * ``CVODE_MEM_NULL`` -- The ``cvode_mem`` pointer is ``NULL``.
     * ``CV_NO_SENS`` -- Sensitivities were not activated.
     * ``CV_NO_QUADSENS`` -- Quadratures depending on the sensitivities were not activated.

   **Notes:**
      When :c:func:`CVodeQuadSensEEtolerances`  is used, before calling :c:func:`CVode`,
      integration of pure quadratures must be initialize and tolerances for pure
      quadratures must be also specified (see :numref:`CVODES.Usage.Purequad`).


.. _CVODES.Usage.FSA.quad.quad_sens_optional_output:

Optional outputs for sensitivity-dependent quadrature integration
-----------------------------------------------------------------

CVODES provides the following functions that can be used to obtain solver
performance information related to quadrature integration.

.. c:function:: int CVodeGetQuadSensNumRhsEvals(void * cvode_mem, long int nrhsQSevals)

   The function :c:func:`CVodeGetQuadSensNumRhsEvals` returns the  number of
   calls made to the user's quadrature right-hand side function.

   **Arguments:**
     * ``cvode_mem`` -- pointer to the CVODES memory block.
     * ``nrhsQSevals`` -- number of calls made to the user's ``rhsQS`` function.

   **Return value:**
     * ``CV_SUCCESS`` -- The optional output value has been successfully set.
     * ``CVODE_MEM_NULL`` -- The ``cvode_mem`` pointer is ``NULL``.
     * ``CV_NO_QUADSENS`` -- Sensitivity-dependent quadrature integration has not been initialized.


.. c:function:: int CVodeGetQuadSensNumErrTestFails(void * cvode_mem, long int nQSetfails)

   The function :c:func:`CVodeGetQuadSensNumErrTestFails` returns the  number of
   local error test failures due to quadrature variables.

   **Arguments:**
     * ``cvode_mem`` -- pointer to the CVODES memory block.
     * ``nQSetfails`` -- number of error test failures due to quadrature variables.

   **Return value:**
     * ``CV_SUCCESS`` -- The optional output value has been successfully set.
     * ``CVODE_MEM_NULL`` -- The ``cvode_mem`` pointer is ``NULL``.
     * ``CV_NO_QUADSENS`` -- Sensitivity-dependent quadrature integration has not been initialized.


.. c:function:: int CVodeGetQuadSensErrWeights(void * cvode_mem, N_Vector * eQSweight)

   The function :c:func:`CVodeGetQuadSensErrWeights` returns the quadrature
   error weights  at the current time.

   **Arguments:**
     * ``cvode_mem`` -- pointer to the CVODES memory block.
     * ``eQSweight`` -- array of quadrature error weight vectors at the current time.

   **Return value:**
     * ``CV_SUCCESS`` -- The optional output value has been successfully set.
     * ``CVODE_MEM_NULL`` -- The ``cvode_mem`` pointer is ``NULL``.
     * ``CV_NO_QUADSENS`` -- Sensitivity-dependent quadrature integration has not been initialized.

   **Notes:**
      .. warning::
         The user must allocate memory for ``eQSweight``.  If quadratures were
         not included in the error control mechanism (through a  call to
         :c:func:`CVodeSetQuadSensErrCon` with ``errconQS = SUNTRUE``), then
         this function does not set the ``eQSweight`` array.


.. c:function:: int CVodeGetQuadSensStats(void * cvode_mem, long int nrhsQSevals, long int nQSetfails)

   The function :c:func:`CVodeGetQuadSensStats` returns the CVODES integrator
   statistics  as a group.

   **Arguments:**
     * ``cvode_mem`` -- pointer to the CVODES memory block.
     * ``nrhsQSevals`` -- number of calls to the user's ``rhsQS`` function.
     * ``nQSetfails`` -- number of error test failures due to quadrature variables.

   **Return value:**
     * ``CV_SUCCESS`` -- the optional output values have been successfully set.
     * ``CVODE_MEM_NULL`` -- the ``cvode_mem`` pointer is ``NULL``.
     * ``CV_NO_QUADSENS`` -- Sensitivity-dependent quadrature integration has not been initialized.


.. _CVODES.Usage.FSA.quad.user_supplied:

User-supplied function for sensitivity-dependent quadrature integration
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

For the integration of sensitivity-dependent quadrature equations, the user must
provide a function that defines the right-hand side of those quadrature
equations. For the sensitivities of quadratures :eq:`CVODES_QUAD` with
integrand :math:`q`, the appropriate right-hand side functions are given by:
:math:`\bar{q}_i = q_y s_i + q_{p_i}`. This user function must be of type
``CVQuadSensRhsFn`` defined as follows:

.. c:type:: int (*CVQuadSensRhsFn)(int Ns, realtype t, N_Vector y, N_Vector *yS, N_Vector yQdot, N_Vector *yQSdot, void *user_data, N_Vector tmp, N_Vector tmpQ)

   This function computes the sensitivity quadrature equation right-hand side
   for a given value of the independent variable :math:`t`` and state vector
   :math:`y`.

   **Arguments:**
      * ``Ns`` -- is the number of sensitivity vectors.
      * ``t`` -- is the current value of the independent variable.
      * ``y`` -- is the current value of the dependent variable vector, :math:`y(t)`.
      * ``ys`` -- is an array of ``Ns`` variables of type ``N_Vector`` containing the dependent sensitivity vectors :math:`s_i`.
      * ``yQdot`` -- is the current value of the quadrature right-hand side, :math:`q`.
      * ``yQSdot``-- array of ``Ns`` vectors to contain the right-hand sides.
      * ``user_data`` -- is the ``user_data`` pointer passed to :c:func:`CVodeSetUserData`.
      * ``tmp1``, ``tmp2`` -- are ``N_Vector`` objects which can be used as temporary storage.

   **Return value:**
      A ``CVQuadSensRhsFn`` should return 0 if successful, a positive value if a
      recoverable error occurred (in which case CVODES will attempt to correct),
      or a negative value if it failed unrecoverably (in which case the
      integration is halted and ``CV_QRHS_FAIL`` is returned).

   **Notes:**
      Allocation of memory for ``rhsvalQS`` is automatically handled within
      CVODES.

      Here ``y`` is of type ``N_Vector`` and ``yS`` is a pointer to an array
      containing ``Ns`` vectors of type ``N_Vector``. It is the user’s
      responsibility to access the vector data consistently (including the use
      of the correct accessor macros from each ``N_Vector`` implementation). For
      the sake of computational efficiency, the vector functions in the two
      ``N_Vector`` implementations provided with CVODES do not perform any
      consistency checks with respect to their ``N_Vector`` arguments.

      There are two situations in which recovery is not possible even if
      ``CVQuadSensRhsFn`` function returns a recoverable error flag. One is when
      this occurs at the very first call to the ``CVQuadSensRhsFn`` (in which
      case CVODES returns ``CV_FIRST_QSRHSFUNC_ERR``). The other is when a
      recoverable error is reported by ``CVQuadSensRhsFn`` after an error test
      failure, while the linear multistep method order is equal to 1 (in which
      case CVODES returns ``CV_UNREC_QSRHSFUNC_ERR``).


.. _CVODES.Usage.FSA.partial:

Note on using partial error control
-----------------------------------

For some problems, when sensitivities are excluded from the error control test,
the behavior of CVODES may appear at first glance to be erroneous. One would
expect that, in such cases, the sensitivity variables would not influence in any
way the step size selection. A comparison of the solver diagnostics reported for
``cvsdenx`` and the second run of the ``cvsfwddenx`` example in
:cite:p:`cvodes_ex` indicates that this may not always be the case.

The short explanation of this behavior is that the step size selection
implemented by the error control mechanism in CVODES is based on the magnitude
of the correction calculated by the nonlinear solver. As mentioned in
:numref:`CVODES.Usage.FSA.user_callable.sensi_malloc`, even with partial error control
selected (in the call to :c:func:`CVodeSetSensErrCon`), the sensitivity
variables are included in the convergence tests of the nonlinear solver.

When using the simultaneous corrector method :numref:`CVODES.Mathematics.FSA`
the nonlinear system that is solved at each step involves both the state and
sensitivity equations. In this case, it is easy to see how the sensitivity
variables may affect the convergence rate of the nonlinear solver and therefore
the step size selection. The case of the staggered corrector approach is more
subtle. After all, in this case (``ism = CV_STAGGERED`` or ``CV_STAGGERED1`` in
the call to :c:func:`CVodeSensInit` :c:func:`CVodeSensInit1`), the sensitivity
variables at a given step are computed only once the solver for the nonlinear
state equations has converged. However, if the nonlinear system corresponding to
the sensitivity equations has convergence problems, CVODES will attempt to
improve the initial guess by reducing the step size in order to provide a better
prediction of the sensitivity variables. Moreover, even if there are no
convergence failures in the solution of the sensitivity system, CVODES may
trigger a call to the linear solver’s setup routine which typically involves
reevaluation of Jacobian information (Jacobian approximation in the case of
CVDENSE and CVBAND, or preconditioner data in the case of the Krylov solvers).
The new Jacobian information will be used by subsequent calls to the nonlinear
solver for the state equations and, in this way, potentially affect the step
size selection.

When using the simultaneous corrector method it is not possible to decide
whether nonlinear solver convergence failures or calls to the linear solver
setup routine have been triggered by convergence problems due to the state or
the sensitivity equations. When using one of the staggered corrector methods
however, these situations can be identified by carefully monitoring the
diagnostic information provided through optional outputs. If there are no
convergence failures in the sensitivity nonlinear solver, and none of the calls
to the linear solver setup routine were made by the sensitivity nonlinear
solver, then the step size selection is not affected by the sensitivity
variables.

Finally, the user must be warned that the effect of appending sensitivity
equations to a given system of ODEs on the step size selection (through the
mechanisms described above) is problem-dependent and can therefore lead to
either an increase or decrease of the total number of steps that CVODES takes to
complete the simulation. At first glance, one would expect that the impact of
the sensitivity variables, if any, would be in the direction of increasing the
step size and therefore reducing the total number of steps. The argument for
this is that the presence of the sensitivity variables in the convergence test
of the nonlinear solver can only lead to additional iterations (and therefore a
smaller final iteration error), or to additional calls to the linear solver
setup routine (and therefore more up-to-date Jacobian information), both of
which will lead to larger steps being taken by CVODES. However, this is true
only locally. Overall, a larger integration step taken at a given time may lead
to step size reductions at later times, due to either nonlinear solver
convergence failures or error test failures.
