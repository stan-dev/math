.. ----------------------------------------------------------------
   SUNDIALS Copyright Start
   Copyright (c) 2002-2022, Lawrence Livermore National Security
   and Southern Methodist University.
   All rights reserved.

   See the top-level LICENSE and NOTICE files for details.

   SPDX-License-Identifier: BSD-3-Clause
   SUNDIALS Copyright End
   ----------------------------------------------------------------

.. _IDAS.Usage.ADJ:

Using IDAS for Adjoint Sensitivity Analysis
===========================================

This chapter describes the use of IDAS to compute sensitivities of derived
functions using adjoint sensitivity analysis. As mentioned before, the adjoint
sensitivity module of IDAS provides the infrastructure for integrating backward
in time any system of DAEs that depends on the solution of the original IVP, by
providing various interfaces to the main IDAS integrator, as well as several
supporting user-callable functions. For this reason, in the following sections
we refer to the *backward problem* and not to the *adjoint problem* when
discussing details relevant to the DAEs that are integrated backward in time.
The backward problem can be the adjoint problem :eq:`IDAS_adj_eqns` or
:eq:`IDAS_adj1_eqns`, and can be augmented with some quadrature differential
equations.

IDAS uses various constants for both input and output. These are defined as
needed in this chapter, but for convenience are also listed separately in
Appendix :numref:`IDAS.Constants`.

We begin with a brief overview, in the form of a skeleton user program.
Following that are detailed descriptions of the interface to the various
user-callable functions and of the user-supplied functions that were not already
described in :numref:`IDAS.Usage.SIM`.

.. _IDAS.Usage.ADJ.skeleton_adj:

A skeleton of the user’s main program
-------------------------------------

The following is a skeleton of the user’s main program as an application of
IDAS. The user program is to have these steps in the order indicated, unless
otherwise noted. For the sake of brevity, we defer many of the details to the
later sections. As in :numref:`IDAS.Usage.SIM.skeleton_sim`, most steps are
independent of the ``N_Vector``, ``SUNMatrix``, ``SUNLinearSolver``, and
``SUNNonlinearSolver`` implementations used. For the steps that are not, refer
to Chapters :numref:`NVectors`, :numref:`SUNMatrix`, :numref:`SUNLinSol`, and
:numref:`SUNNonlinSol` for the specific name of the function to be called or
macro to be referenced.

Steps that are changed from the skeleton programs presented in
:numref:`IDAS.Usage.SIM.skeleton_sim`, :numref:`IDAS.Usage.FSA.skeleton_sim`,
and :numref:`IDAS.Usage.FSA.quad`, are bolded.

#. Initialize parallel or multi-threaded environment

#. Create the SUNDIALS context object

**Forward Problem**

#. Set initial conditions for the forward problem

#. Create matrix object for the forward problem

#. Create linear solver object for the forward problem

#. Create nonlinear solver module for the forward problem

#. Create IDAS object for the forward problem

#. Initialize IDAS solver for the forward problem

#. Specify integration tolerances for forward problem

#. Attach linear solver module for the forward problem

#. Set linear solver optional inputs for the forward problem

#. Attach nonlinear solver module for the forward problem

#. Set nonlinear solver optional inputs for the forward problem

#. Initialize quadrature problem or problems for forward problems, using :c:func:`IDAQuadInit` and/or :c:func:`IDAQuadSensInit`.

#. Initialize forward sensitivity problem

#. Specify rootfinding

#. Set optional inputs for the forward problem

#. **Allocate space for the adjoint computation**

   Call :c:func:`IDAAdjInit` to allocate memory for the combined
   forward-backward problem. This call requires ``Nd``, the number of steps
   between two consecutive checkpoints. :c:func:`IDAAdjInit` also specifies the
   type of interpolation used (see :numref:`IDAS.Mathematics.ASA.Checkpointing`).

#. **Integrate forward problem**

   Call :c:func:`IDASolveF`, a wrapper for the IDAS main integration function
   :c:func:`IDASolve`, either in ``IDA_NORMAL`` mode to the time ``tout`` or in
   ``IDA_ONE_STEP`` mode inside a loop (if intermediate solutions of the forward
   problem are desired (see :numref:`IDAS.Usage.ADJ.user_callable.idasolvef`).
   The final value of ``tret`` is then the maximum allowable value for the
   endpoint :math:`T` of the backward problem.

**Backward Problem(s)**

.. _IDAS.Usage.ADJ.skeleton_adj.back_start:

18. **Create vectors of endpoint values for the backward problem**

    Create the vectors ``yB0`` and ``ypB0`` at the endpoint time
    ``tB0`` :math:`= T` at which the backward problem starts.

#.  **Create the backward problem**

    Call :c:func:`IDACreateB`, a wrapper for :c:func:`IDACreate`, to create the
    IDAS memory block for the new backward problem. Unlike :c:func:`IDACreate`,
    the function :c:func:`IDACreateB` does not return a pointer to the newly
    created memory block (see :numref:`IDAS.Usage.ADJ.user_callable.idainitb`).
    Instead, this pointer is attached to the internal adjoint memory block
    (created by :c:func:`IDAAdjInit`) and returns an identifier called ``which``
    that the user must later specify in any actions on the newly created backward
    problem.

#.  **Allocate memory for the backward problem**

    Call :c:func:`IDAInitB` (or :c:func:`IDAInitBS`, when the backward problem
    depends on the forward sensitivities). The two functions are actually
    wrappers for :c:func:`IDAInit` and allocate internal memory, specify problem
    data, and initialize IDAS at ``tB0`` for the backward problem (see
    :numref:`IDAS.Usage.ADJ.user_callable.idainitb`).

#.  **Specify integration tolerances for backward problem**

    Call :c:func:`IDASStolerancesB` or :c:func:`IDASVtolerancesB` to specify a
    scalar relative tolerance and scalar absolute tolerance, or a scalar relative
    tolerance and a vector of absolute tolerances, respectively. The functions
    are wrappers for :c:func:`IDASStolerances` and :c:func:`IDASVtolerances` but
    they require an extra argument ``which``, the identifier of the backward
    problem returned by :c:func:`IDACreateB`. See
    :numref:`IDAS.Usage.ADJ.user_callable.idatolerances_b` for more information.

#.  **Set optional inputs for the backward problem**

    Call ``IDASet*B`` functions to change from their default values any optional
    inputs that control the behavior of IDAS. Unlike their counterparts for the
    forward problem, these functions take an extra argument ``which``, the
    identifier of the backward problem returned by :c:func:`IDACreateB` (see
    :numref:`IDAS.Usage.ADJ.user_callable.optional_input_b`).


.. _IDAS.Usage.ADJ.skeleton_adj.matrixB:

23. **Create matrix object for the backward problem**

    If a nonlinear solver requiring a linear solve will be used (e.g., the the
    default Newton iteration) and the linear solver will be a direct linear
    solver, then a template Jacobian matrix must be created by calling the
    appropriate constructor function defined by the particular ``SUNMatrix``
    implementation.

    .. note::

       The dense, banded, and sparse matrix objects are usable only in a serial
       or threaded environment.

       It is not required to use the same matrix type for both the forward and
       the backward problems.


.. _IDAS.Usage.ADJ.skeleton_adj.lin_solverB:

24. **Create linear solver object for the backward problem**

    If a nonlinear solver requiring a linear solver is chosen (e.g., the default
    Newton iteration), then the desired linear solver object for the backward
    problem must be created by calling the appropriate constructor function
    defined by the particular ``SUNLinearSolver`` implementation.

    .. note::

       It is not required to use the same linear solver module for both the
       forward and the backward problems; for example, the forward problem could
       be solved with the ``SUNLINSOL_BAND`` linear solver module and the
       backward problem with ``SUNLINSOL_SPGMR`` linear solver module.

#.  **Set linear solver interface optional inputs for the backward problem**

    Call ``IDASet*B`` functions to change optional inputs specific to the linear
    solver interface. See :numref:`IDAS.Usage.ADJ.user_callable.optional_input_b`
    for details.

#.  **Attach linear solver module for the backward problem**

    If a nonlinear solver requiring a linear solver is chosen for the backward
    problem (e.g., the default Newton iteration), then initialize the IDALS
    linear solver interface by attaching the linear solver object (and matrix
    object, if applicable) with :c:func:`IDASetLinearSolverB` (for additional
    details see :numref:`IDAS.Usage.ADJ.user_callable.lin_solv_b`).

#.  **Create nonlinear solver object for the backward problem** (*optional*)

    If using a non-default nonlinear solver for the backward problem, then create
    the desired nonlinear solver object by calling the appropriate constructor
    function defined by the particular ``SUNNonlinearSolver`` implementation
    e.g., ``NLSB = SUNNonlinSol_***(...);`` where ``***`` is the name of the
    nonlinear solver (see Chapter :numref:`SUNNonlinSol` for details).

#.  **Attach nonlinear solver module for the backward problem** (*optional*)

    If using a non-default nonlinear solver for the backward problem, then
    initialize the nonlinear solver interface by attaching the nonlinear solver
    object by calling :c:func:`IDASetNonlinearSolverB`.


.. _IDAS.Usage.ADJ.skeleton_adj.quadB:

29. **Initialize quadrature calculation**

    If additional quadrature equations must be evaluated, call
    :c:func:`IDAQuadInitB` or :c:func:`IDAQuadInitBS` (if quadrature depends
    also on the forward sensitivities) as shown in
    :numref:`IDAS.Usage.ADJ.user_callable.backquad.idaquadinitb`. These
    functions are wrappers around :c:func:`IDAQuadInit` and can be used to
    initialize and allocate memory for quadrature integration. Optionally, call
    ``IDASetQuad*B`` functions to change from their default values optional
    inputs that control the integration of quadratures during the backward
    phase.

#.  **Integrate backward problem**

    Call :c:func:`IDASolveB`, a second wrapper around the IDAS main integration
    function :c:func:`IDASolve`, to integrate the backward problem from ``tB0``.
    This function can be called either in ``IDA_NORMAL`` or ``IDA_ONE_STEP``
    mode. Typically, :c:func:`IDASolveB` will be called in ``IDA_NORMAL`` mode
    with an end time equal to the initial time :math:`t_0` of the forward
    problem.


.. _IDAS.Usage.ADJ.skeleton_adj.back_end:

31. **Extract quadrature variables**

    If applicable, call :c:func:`IDAGetQuadB`, a wrapper around
    :c:func:`IDAGetQuad`, to extract the values of the quadrature variables at
    the time returned by the last call to :c:func:`IDASolveB`.

#.  **Deallocate memory**

    Upon completion of the backward integration, call all necessary deallocation
    functions. These include appropriate destructors for the vectors ``y`` and
    ``yB``, a call to :c:func:`IDAFree` to free the IDAS memory block for the
    forward problem. If one or more additional adjoint sensitivity analyses are
    to be done for this problem, a call to :c:func:`IDAAdjFree` (see
    :numref:`IDAS.Usage.ADJ.user_callable.idaadjinit`) may be made to free and
    deallocate the memory allocated for the backward problems, followed by a call
    to :c:func:`IDAAdjInit`.

#.  Finalize MPI, if used

The above user interface to the adjoint sensitivity module in IDAS was motivated
by the desire to keep it as close as possible in look and feel to the one for
DAE IVP integration. Note that if steps
(:ref:`18<IDAS.Usage.ADJ.skeleton_adj.back_start>`) -
(:ref:`31<IDAS.Usage.ADJ.skeleton_adj.back_end>`) are not present, a program with
the above structure will have the same functionality as one described in
:numref:`IDAS.Usage.SIM.skeleton_sim` for integration of DAEs, albeit with some
overhead due to the checkpointing scheme.

If there are multiple backward problems associated with the same forward
problem, repeat steps (:ref:`18<IDAS.Usage.ADJ.skeleton_adj.back_start>`) -
(:ref:`31<IDAS.Usage.ADJ.skeleton_adj.back_end>`) above for each successive
backward problem. In the process, If there are multiple backward problems
associated with the same forward each call to :c:func:`IDACreateB` creates a new
value of the identifier ``which``.


.. _IDAS.Usage.ADJ.user_callable:

User-callable functions for adjoint sensitivity analysis
--------------------------------------------------------

.. _IDAS.Usage.ADJ.user_callable.idaadjinit:

Adjoint sensitivity allocation and deallocation functions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

After the setup phase for the forward problem, but before the call to
:c:func:`IDASolveF`, memory for the combined forward-backward problem must be
allocated by a call to the function :c:func:`IDAAdjInit`. The form of the call
to this function is

.. c:function:: int IDAAdjInit(void * ida_mem, long int Nd, int interpType)

   The function :c:func:`IDAAdjInit` updates IDAS memory block by allocating
   the internal memory needed for backward integration.  Space is allocated for
   the ``Nd`` :math:`= N_d` interpolation data points, and  a linked list of
   checkpoints is initialized.

   **Arguments:**
     * ``ida_mem`` -- is the pointer to the IDAS memory block returned by a previous call to :c:func:`IDACreate`.
     * ``Nd`` -- is the number of integration steps between two consecutive checkpoints.
     * ``interpType`` -- specifies the type of interpolation used and can be ``IDA_POLYNOMIAL`` or ``IDA_HERMITE`` , indicating variable-degree polynomial and cubic Hermite interpolation, respectively see :numref:`IDAS.Mathematics.ASA.Checkpointing`.

   **Return value:**
     * ``IDA_SUCCESS`` -- :c:func:`IDAAdjInit` was successful.
     * ``IDA_MEM_FAIL`` -- A memory allocation request has failed.
     * ``IDA_MEM_NULL`` -- ``ida_mem`` was ``NULL``.
     * ``IDA_ILL_INPUT`` -- One of the parameters was invalid: ``Nd`` was not positive or ``interpType`` is not one of the ``IDA_POLYNOMIAL`` or ``IDA_HERMITE``.

   **Notes:**

   The user must set ``Nd`` so that all data needed for interpolation of the
   forward problem solution between two checkpoints fits in memory.
   :c:func:`IDAAdjInit` attempts to allocate space for :math:`(2 N_d+3)`
   variables of type ``N_Vector``.

   If an error occurred, :c:func:`IDAAdjInit` also sends a message to the error
   handler function.


.. c:function:: int IDAAdjReInit(void * ida_mem)

   The function :c:func:`IDAAdjReInit` reinitializes the IDAS memory  block for
   ASA, assuming that the number of steps between check  points and the type of
   interpolation remain unchanged.

   **Arguments:**
     * ``ida_mem`` -- is the pointer to the IDAS memory block returned by a previous call to :c:func:`IDACreate`.

   **Return value:**
     * ``IDA_SUCCESS`` -- :c:func:`IDAAdjReInit` was successful.
     * ``IDA_MEM_NULL`` -- ``ida_mem`` was ``NULL``.
     * ``IDA_NO_ADJ`` -- The function :c:func:`IDAAdjInit` was not previously called.

   **Notes:**

   The list of check points (and associated memory) is deleted.

   The list of backward problems is kept. However, new backward problems can be
   added to this list by calling :c:func:`IDACreateB`. If a new list of backward
   problems is also needed, then free the adjoint memory (by calling
   :c:func:`IDAAdjFree`) and reinitialize ASA with :c:func:`IDAAdjInit`.

   The IDAS memory for the forward and backward problems can be reinitialized
   separately by calling :c:func:`IDAReInit` and :c:func:`IDAReInitB`,
   respectively.


.. c:function:: void IDAAdjFree(void * ida_mem)

   The function :c:func:`IDAAdjFree` frees the memory related to backward integration
   allocated by a previous call to :c:func:`IDAAdjInit`.


   **Arguments:**
      The only argument is the IDAS memory block pointer returned by a previous
      call to :c:func:`IDACreate`.

   **Return value:**
      The function :c:func:`IDAAdjFree` has no return value.

   **Notes:**

   This function frees all memory allocated by :c:func:`IDAAdjInit`. This
   includes workspace memory, the linked list of checkpoints, memory for the
   interpolation data, as well as the IDAS memory for the backward integration
   phase.

   Unless one or more further calls to :c:func:`IDAAdjInit` are to be made,
   :c:func:`IDAAdjFree` should not be called by the user, as it is invoked
   automatically by :c:func:`IDAFree`.


Adjoint sensitivity optional input
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

At any time during the integration of the forward problem, the user can disable
the checkpointing of the forward sensitivities by calling the following
function:

.. c:function:: int IDAAdjSetNoSensi(void * ida_mem)

   The function :c:func:`IDAAdjSetNoSensi` instructs :c:func:`IDASolveF` not  to
   save checkpointing data for forward sensitivities any more.

   **Arguments:**
     * ``ida_mem`` -- pointer to the IDAS memory block.

   **Return value:**
     * ``IDA_SUCCESS`` -- The call to :c:func:`IDACreateB` was successful.
     * ``IDA_MEM_NULL`` -- The ``ida_mem`` was ``NULL``.
     * ``IDA_NO_ADJ`` -- The function :c:func:`IDAAdjInit` has not been previously called.


.. _IDAS.Usage.ADJ.user_callable.idasolvef:

Forward integration function
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The function :c:func:`IDASolveF` is very similar to the IDAS function
:c:func:`IDASolve` in that it integrates the solution of the forward problem and
returns the solution :math:`(y,\dot{y})`. At the same time, however,
:c:func:`IDASolveF` stores checkpoint data every ``Nd`` integration steps.
:c:func:`IDASolveF` can be called repeatedly by the user. Note that
:c:func:`IDASolveF` is used only for the forward integration pass within an
Adjoint Sensitivity Analysis. It is not for use in Forward Sensitivity Analysis;
for that, see :numref:`IDAS.Usage.FSA`. The call to this function
has the form

.. c:function:: int IDASolveF(void * ida_mem, realtype tout, realtype * tret, N_Vector yret, N_Vector ypret, int itask, int * ncheck)

   The function :c:func:`IDASolveF` integrates the forward problem over an
   interval in :math:`t`  and saves checkpointing data.

   **Arguments:**
     * ``ida_mem`` -- pointer to the IDAS memory block.
     * ``tout`` -- the next time at which a computed solution is desired.
     * ``tret`` -- the time reached by the solver output.
     * ``yret`` -- the computed solution vector :math:`y`.
     * ``ypret`` -- the computed solution vector :math:`\dot{y}`.
     * ``itask`` --  a flag indicating the job of the solver for the next step. The ``IDA_NORMAL`` task is to have the solver take internal steps until it has reached or just passed the user-specified ``tout`` parameter. The solver then interpolates in order to return an approximate value of :math:`y(\texttt{tout})` and :math:`\dot{y}(\texttt{tout})`. The ``IDA_ONE_STEP`` option tells the solver to take just one internal step and return the solution at the point reached by that step.
     * ``ncheck`` -- the number of internal checkpoints stored so far.

   **Return value:**

   On return, :c:func:`IDASolveF` returns vectors ``yret``, ``ypret`` and a
   corresponding independent variable value ``t = tret``, such that ``yret`` is
   the computed value of :math:`y(t)` and ``ypret`` the value of
   :math:`\dot{y}(t)`. Additionally, it returns in ``ncheck`` the number of
   internal checkpoints saved; the total number of checkpoint intervals is
   ``ncheck+1``. The return value flag (of type ``int``) will be one of the
   following. For more details see the documentation for :c:func:`IDASolve`.

     * ``IDA_SUCCESS`` -- :c:func:`IDASolveF` succeeded.
     * ``IDA_TSTOP_RETURN`` -- :c:func:`IDASolveF` succeeded by reaching the optional stopping point.
     * ``IDA_ROOT_RETURN`` -- :c:func:`IDASolveF` succeeded and found one or more roots. In this case, ``tret`` is the location of the root. If ``nrtfn`` :math:`>1` , call :c:func:`IDAGetRootInfo` to see which :math:`g_i` were found to have a root.
     * ``IDA_NO_MALLOC`` -- The function :c:func:`IDAInit` has not been previously called.
     * ``IDA_ILL_INPUT`` -- One of the inputs to :c:func:`IDASolveF` is illegal.
     * ``IDA_TOO_MUCH_WORK`` -- The solver took ``mxstep`` internal steps but could not reach ``tout``.
     * ``IDA_TOO_MUCH_ACC`` -- The solver could not satisfy the accuracy demanded by the user for some internal step.
     * ``IDA_ERR_FAILURE`` -- Error test failures occurred too many times during one internal time step or occurred with :math:`|h| = h_{min}`.
     * ``IDA_CONV_FAILURE`` -- Convergence test failures occurred too many times during one internal time step or occurred with :math:`|h| = h_{min}`.
     * ``IDA_LSETUP_FAIL`` -- The linear solver's setup function failed in an unrecoverable manner.
     * ``IDA_LSOLVE_FAIL`` -- The linear solver's solve function failed in an unrecoverable manner.
     * ``IDA_NO_ADJ`` -- The function :c:func:`IDAAdjInit` has not been previously called.
     * ``IDA_MEM_FAIL`` -- A memory allocation request has failed in an attempt to allocate space for a new checkpoint.

   **Notes:**

   All failure return values are negative and therefore a test ``flag``:math:`<
   0` will trap all :c:func:`IDASolveF` failures.

   At this time, :c:func:`IDASolveF` stores checkpoint information in memory
   only.  Future versions will provide for a safeguard option of dumping
   checkpoint data into a temporary file as needed. The data stored at each
   checkpoint is basically a snapshot of the IDAS internal memory block and
   contains enough information to restart the integration from that time and to
   proceed with the same step size and method order sequence as during the
   forward integration.

   In addition, :c:func:`IDASolveF` also stores interpolation data between
   consecutive checkpoints so that, at the end of this first forward integration
   phase, interpolation information is already available from the last
   checkpoint forward. In particular, if no checkpoints were necessary, there is
   no need for the second forward integration phase.

   .. warning::

      It is illegal to change the integration tolerances between consecutive
      calls  to :c:func:`IDASolveF`, as this information is not captured in
      the checkpoint data.


.. _IDAS.Usage.ADJ.user_callable.idainitb:

Backward problem initialization functions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The functions :c:func:`IDACreateB` and :c:func:`IDAInitB` (or :c:func:`IDAInitBS`) must be called
in the order listed. They instantiate an IDAS solver object, provide problem and
solution specifications, and allocate internal memory for the backward problem.

.. c:function:: int IDACreateB(void * ida_mem, int * which)

   The function :c:func:`IDACreateB` instantiates an IDAS solver object for the
   backward problem.

   **Arguments:**
     * ``ida_mem`` -- pointer to the IDAS memory block returned by :c:func:`IDACreate`.
     * ``which`` -- contains the identifier assigned by IDAS for the newly created backward problem. Any call to ``IDA*B`` functions requires such an identifier.

   **Return value:**
     * ``IDA_SUCCESS`` -- The call to :c:func:`IDACreateB` was successful.
     * ``IDA_MEM_NULL`` -- The ``ida_mem`` was ``NULL``.
     * ``IDA_NO_ADJ`` -- The function :c:func:`IDAAdjInit` has not been previously called.
     * ``IDA_MEM_FAIL`` -- A memory allocation request has failed.


There are two initialization functions for the backward problem – one for the
case when the backward problem does not depend on the forward sensitivities, and
one for the case when it does. These two functions are described next.

The function :c:func:`IDAInitB` initializes the backward problem when it does
not depend on the forward sensitivities. It is essentially wrapper for IDAInit
with some particularization for backward integration, as described below.


.. c:function:: int IDAInitB(void * ida_mem, int which, IDAResFnB resB, realtype tB0, N_Vector yB0, N_Vector ypB0)

   The function :c:func:`IDAInitB` provides problem specification, allocates
   internal memory,  and initializes the backward problem.

   **Arguments:**
     * ``ida_mem`` -- pointer to the IDAS memory block returned by :c:func:`IDACreate`.
     * ``which`` -- represents the identifier of the backward problem.
     * ``resB`` -- is the C function which computes :math:`fB` , the residual of the backward DAE problem. This function has the form ``resB(t, y, yp, yB, ypB, resvalB, user_dataB)`` for full details see :numref:`IDAS.Usage.ADJ.user_supplied.DAEres_b`.
     * ``tB0`` -- specifies the endpoint :math:`T` where final conditions are provided for the backward problem, normally equal to the endpoint of the forward integration.
     * ``yB0`` -- is the initial value at :math:`t =` ``tB0`` of the backward solution.
     * ``ypB0`` -- is the initial derivative value at :math:`t =` ``tB0`` of the backward solution.

   **Return value:**
     * ``IDA_SUCCESS`` -- The call to :c:func:`IDAInitB` was successful.
     * ``IDA_NO_MALLOC`` -- The function :c:func:`IDAInit` has not been previously called.
     * ``IDA_MEM_NULL`` -- The ``ida_mem`` was ``NULL``.
     * ``IDA_NO_ADJ`` -- The function :c:func:`IDAAdjInit` has not been previously called.
     * ``IDA_BAD_TB0`` -- The final time ``tB0`` was outside the interval over which the forward problem was solved.
     * ``IDA_ILL_INPUT`` -- The parameter ``which`` represented an invalid identifier, or one of ``yB0`` , ``ypB0`` , ``resB`` was ``NULL``.

   **Notes:**
      The memory allocated by :c:func:`IDAInitB` is deallocated by the function  :c:func:`IDAAdjFree`.


For the case when backward problem also depends on the forward sensitivities,
user must call :c:func:`IDAInitBS` instead of :c:func:`IDAInitB`. Only the third
argument of each function differs between these functions.


.. c:function:: int IDAInitBS(void * ida_mem, int which, IDAResFnBS resBS, realtype tB0, N_Vector yB0, N_Vector ypB0)

   The function :c:func:`IDAInitBS` provides problem specification, allocates
   internal memory,  and initializes the backward problem.

   **Arguments:**
     * ``ida_mem`` -- pointer to the IDAS memory block returned by :c:func:`IDACreate`.
     * ``which`` -- represents the identifier of the backward problem.
     * ``resBS`` -- is the C function which computes :math:`fB` , the residual or the backward DAE problem. This function has the form ``resBS(t, y, yp, yS, ypS, yB, ypB, resvalB, user_dataB)`` for full details see :numref:`IDAS.Usage.ADJ.DAEres_bs`.
     * ``tB0`` -- specifies the endpoint :math:`T` where final conditions are provided for the backward problem.
     * ``yB0`` -- is the initial value at :math:`t =` ``tB0`` of the backward solution.
     * ``ypB0`` -- is the initial derivative value at :math:`t =` ``tB0`` of the backward solution.

   **Return value:**
     * ``IDA_SUCCESS`` -- The call to :c:func:`IDAInitB` was successful.
     * ``IDA_NO_MALLOC`` -- The function :c:func:`IDAInit` has not been previously called.
     * ``IDA_MEM_NULL`` -- The ``ida_mem`` was ``NULL``.
     * ``IDA_NO_ADJ`` -- The function :c:func:`IDAAdjInit` has not been previously called.
     * ``IDA_BAD_TB0`` -- The final time ``tB0`` was outside the interval over which the forward problem was solved.
     * ``IDA_ILL_INPUT`` -- The parameter ``which`` represented an invalid identifier, or one of ``yB0`` , ``ypB0`` , ``resB`` was ``NULL`` , or sensitivities were not active during the forward integration.

   **Notes:**
      The memory allocated by :c:func:`IDAInitBS` is deallocated by the function
      :c:func:`IDAAdjFree`.


The function :c:func:`IDAReInitB` reinitializes idas for the solution of a
series of backward problems, each identified by a value of the parameter which.
:c:func:`IDAReInitB` is essentially a wrapper for :c:func:`IDAReInit`, and so
all details given for :c:func:`IDAReInit` apply here. Also, :c:func:`IDAReInitB`
can be called to reinitialize a backward problem even if it has been initialized
with the sensitivity-dependent version :c:func:`IDAInitBS`. Before calling
:c:func:`IDAReInitB` for a new backward problem, call any desired solution
extraction functions ``IDAGet**`` associated with the previous backward problem.
The call to the :c:func:`IDAReInitB` function has the form


.. c:function:: int IDAReInitB(void * ida_mem, int which, realtype tB0, N_Vector yB0, N_Vector ypB0)

   The function :c:func:`IDAReInitB` reinitializes an IDAS backward problem.

   **Arguments:**
     * ``ida_mem`` -- pointer to IDAS memory block returned by :c:func:`IDACreate`.
     * ``which`` -- represents the identifier of the backward problem.
     * ``tB0`` -- specifies the endpoint :math:`T` where final conditions are provided for the backward problem.
     * ``yB0`` -- is the initial value at :math:`t =` ``tB0`` of the backward solution.
     * ``ypB0`` -- is the initial derivative value at :math:`t =` ``tB0`` of the backward solution.

   **Return value:**
     * ``IDA_SUCCESS`` -- The call to :c:func:`IDAReInitB` was successful.
     * ``IDA_NO_MALLOC`` -- The function :c:func:`IDAInit` has not been previously called.
     * ``IDA_MEM_NULL`` -- The ``ida_mem`` memory block pointer was ``NULL``.
     * ``IDA_NO_ADJ`` -- The function :c:func:`IDAAdjInit` has not been previously called.
     * ``IDA_BAD_TB0`` -- The final time ``tB0`` is outside the interval over which the forward problem was solved.
     * ``IDA_ILL_INPUT`` -- The parameter ``which`` represented an invalid identifier, or one of ``yB0`` , ``ypB0`` was ``NULL``.


.. _IDAS.Usage.ADJ.user_callable.idatolerances_b:

Tolerance specification functions for backward problem
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

One of the following two functions must be called to specify the integration
tolerances for the backward problem. Note that this call must be made after the
call to :c:func:`IDAInitB` or :c:func:`IDAInitBS`.

.. c:function:: int IDASStolerancesB(void * ida_mem, int which, realtype reltolB, realtype abstolB)

   The function :c:func:`IDASStolerancesB` specifies scalar relative and
   absolute  tolerances.

   **Arguments:**
     * ``ida_mem`` -- pointer to the IDAS memory block returned by :c:func:`IDACreate`.
     * ``which`` -- represents the identifier of the backward problem.
     * ``reltolB`` -- is the scalar relative error tolerance.
     * ``abstolB`` -- is the scalar absolute error tolerance.

   **Return value:**
     * ``IDA_SUCCESS`` -- The call to :c:func:`IDASStolerancesB` was successful.
     * ``IDA_MEM_NULL`` -- The IDAS memory block was not initialized through a previous call to :c:func:`IDACreate`.
     * ``IDA_NO_MALLOC`` -- The allocation function :c:func:`IDAInit` has not been called.
     * ``IDA_NO_ADJ`` -- The function :c:func:`IDAAdjInit` has not been previously called.
     * ``IDA_ILL_INPUT`` -- One of the input tolerances was negative.


.. c:function:: int IDASVtolerancesB(void * ida_mem, int which, realtype reltolB, N_Vector abstolB)

   The function :c:func:`IDASVtolerancesB` specifies scalar relative tolerance
   and  vector absolute tolerances.

   **Arguments:**
     * ``ida_mem`` -- pointer to the IDAS memory block returned by :c:func:`IDACreate`.
     * ``which`` -- represents the identifier of the backward problem.
     * ``reltolB`` -- is the scalar relative error tolerance.
     * ``abstolB`` -- is the vector of absolute error tolerances.

   **Return value:**
     * ``IDA_SUCCESS`` -- The call to :c:func:`IDASVtolerancesB` was successful.
     * ``IDA_MEM_NULL`` -- The IDAS memory block was not initialized through a previous call to :c:func:`IDACreate`.
     * ``IDA_NO_MALLOC`` -- The allocation function :c:func:`IDAInit` has not been called.
     * ``IDA_NO_ADJ`` -- The function :c:func:`IDAAdjInit` has not been previously called.
     * ``IDA_ILL_INPUT`` -- The relative error tolerance was negative or the absolute tolerance had a negative component.

   **Notes:**
      This choice of tolerances is important when the absolute error tolerance
      needs to  be different for each component of the DAE state vector
      :math:`y`.


.. _IDAS.Usage.ADJ.user_callable.lin_solv_b:

Linear solver initialization functions for backward problem
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

All IDAS linear solver modules available for forward problems are available for
the backward problem. They should be created as for the forward problem then
attached to the memory structure for the backward problem using the following
function.

.. c:function:: int IDASetLinearSolverB(void * ida_mem, int which, SUNLinearSolver LS, SUNMatrix A)

   The function :c:func:`IDASetLinearSolverB` attaches a generic
   ``SUNLinearSolver`` object ``LS`` and corresponding template  Jacobian
   ``SUNMatrix`` object ``A`` (if applicable) to IDAS,  initializing the IDALS
   linear solver interface for solution of  the backward problem.

   **Arguments:**
     * ``ida_mem`` -- pointer to the IDAS memory block.
     * ``which`` -- represents the identifier of the backward problem returned by :c:func:`IDACreateB`.
     * ``LS`` -- SUNLinearSolver object to use for solving linear systems for the backward problem.
     * ``A`` -- SUNMatrix object for used as a template for the Jacobian for the backward problem or ``NULL`` if not applicable.

   **Return value:**
     * ``IDALS_SUCCESS`` -- The IDALS initialization was successful.
     * ``IDALS_MEM_NULL`` -- The ``ida_mem`` pointer is ``NULL``.
     * ``IDALS_ILL_INPUT`` -- The parameter ``which`` represented an invalid identifier.
     * ``IDALS_MEM_FAIL`` -- A memory allocation request failed.
     * ``IDALS_NO_ADJ`` -- The function :c:func:`IDAAdjInit` has not been previously called.

   **Notes:**

   If ``LS`` is a matrix-based linear solver, then the template Jacobian matrix
   ``A`` will be used in the solve process, so if additional storage is required
   within the ``SUNMatrix`` object (e.g. for factorization of a banded matrix),
   ensure that the input object is allocated with sufficient size (see the
   documentation of the particular ``SUNMatrix`` type in Chapter
   :numref:`SUNMatrix` for further information).

   The previous routines ``IDADlsSetLinearSolverB`` and
   ``IDASpilsSetLinearSolverB`` are now deprecated.


.. _IDAS.Usage.ADJ.user_callable.nonlin_solv_init_b:

Nonlinear solver initialization functions for backward problem
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

As with the forward problem IDAS uses the ``SUNNonlinearSolver`` implementation
of Newton’s method defined by the ``SUNNONLINSOL_NEWTON`` module (see
:numref:`SUNNonlinSol.Newton`) by default.

To specify a different nonlinear solver in IDAS for the backward problem, the
user’s program must create a ``SUNNonlinearSolver`` object by calling the
appropriate constructor routine. The user must then attach the
``SUNNonlinearSolver`` object to IDAS by calling
:c:func:`IDASetNonlinearSolverB`, as documented below.

When changing the nonlinear solver in IDAS, :c:func:`IDASetNonlinearSolverB`
must be called after :c:func:`IDAInitB`. If any calls to :c:func:`IDASolveB`
have been made, then IDAS will need to be reinitialized by calling
:c:func:`IDAReInitB` to ensure that the nonlinear solver is initialized
correctly before any subsequent calls to :c:func:`IDASolveB`.


.. c:function:: int IDASetNonlinearSolverB(void * ida_mem, int which, SUNNonlinearSolver NLS)

   The function :c:func:`IDASetNonLinearSolverB` attaches a
   ``SUNNonlinearSolver``  object (``NLS``) to IDAS for the solution of the
   backward problem.

   **Arguments:**
     * ``ida_mem`` -- pointer to the IDAS memory block.
     * ``which`` -- represents the identifier of the backward problem returned by :c:func:`IDACreateB`.
     * ``NLS`` -- SUNNonlinearSolver object to use for solving nonlinear systems for the backward problem.

   **Return value:**
     * ``IDA_SUCCESS`` -- The nonlinear solver was successfully attached.
     * ``IDA_MEM_NULL`` -- The ``ida_mem`` pointer is ``NULL``.
     * ``IDALS_NO_ADJ`` -- The function ``IDAAdjInit`` has not been previously called.
     * ``IDA_ILL_INPUT`` -- The parameter ``which`` represented an invalid identifier or the SUNNonlinearSolver object is ``NULL`` , does not implement the required nonlinear solver operations, is not of the correct type, or the residual function, convergence test function, or maximum number of nonlinear iterations could not be set.


.. _sss:idacalcicB:

Initial condition calculation functions for backward problem
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

IDAS provides support for calculation of consistent initial conditions for
certain backward index-one problems of semi-implicit form through the functions
:c:func:`IDACalcICB` and :c:func:`IDACalcICBS`. Calling them is optional. It is
only necessary when the initial conditions do not satisfy the adjoint system.

The above functions provide the same functionality for backward problems as
:c:func:`IDACalcIC` with parameter ``icopt`` = ``IDA_YA_YDP_INIT`` provides for
forward problems: compute the algebraic components of :math:`yB` and
differential components of :math:`\dot{y}B`, given the differential components
of :math:`yB`. They require that the :c:func:`IDASetIdB` was previously called
to specify the differential and algebraic components.

Both functions require forward solutions at the final time ``tB0``.
:c:func:`IDACalcICBS` also needs forward sensitivities at the final time
``tB0``.

.. c:function:: int IDACalcICB(void * ida_mem, int which, realtype tBout1, N_Vector yfin, N_Vector ypfin)

   The function :c:func:`IDACalcICB` corrects the initial values ``yB0`` and
   ``ypB0`` at  time ``tB0`` for the backward problem.

   **Arguments:**
     * ``ida_mem`` -- pointer to the IDAS memory block.
     * ``which`` -- is the identifier of the backward problem.
     * ``tBout1`` -- is the first value of :math:`t` at which a solution will be requested from :c:func:`IDASolveB`. This value is needed here only to determine the direction of integration and rough scale in the independent variable :math:`t`.
     * ``yfin`` -- the forward solution at the final time ``tB0``.
     * ``ypfin`` -- the forward solution derivative at the final time ``tB0``.

   **Return value:**
     * ``IDA_NO_ADJ`` -- :c:func:`IDAAdjInit` has not been previously called.
     * ``IDA_ILL_INPUT`` -- Parameter ``which`` represented an invalid identifier.

   **Notes:**

   All failure return values are negative and therefore a test ``flag`` :math:`<
   0` will trap all :c:func:`IDACalcICB` failures.  Note that
   :c:func:`IDACalcICB` will correct the values of :math:`yB(tB_0)` and
   :math:`\dot{y}B(tB_0)` which were specified in the previous call to
   :c:func:`IDAInitB` or :c:func:`IDAReInitB`. To obtain the corrected values,
   call :c:func:`IDAGetconsistentICB` (see
   :numref:`IDAS.Usage.ADJ.user_callable.optional_ouput_b.iccalcB`).

   :c:func:`IDACalcICB` will correct the values of :math:`yB(tB_0)` and
   :math:`\dot{y}B(tB_0)` which were specified in the previous call to
   :c:func:`IDAInitB` or :c:func:`IDAReInitB`. To obtain the corrected values,
   :call c:func:`IDAGetConsistentICB` (see
   ::numref:`IDAS.Usage.ADJ.user_callable.optional_output_b`).

In the case where the backward problem also depends on the forward
sensitivities, user must call the following function to correct the initial
conditions:

.. c:function:: int IDACalcICBS(void * ida_mem, int which, realtype tBout1, N_Vector yfin, N_Vector ypfin, N_Vector ySfin, N_Vector ypSfin)

   The function :c:func:`IDACalcICBS` corrects the initial values ``yB0`` and
   ``ypB0`` at  time ``tB0`` for the backward problem.

   **Arguments:**
     * ``ida_mem`` -- pointer to the IDAS memory block.
     * ``which`` -- is the identifier of the backward problem.
     * ``tBout1`` -- is the first value of :math:`t` at which a solution will be requested from :c:func:`IDASolveB` .This value is needed here only to determine the direction of integration and rough scale in the independent variable :math:`t`.
     * ``yfin`` -- the forward solution at the final time ``tB0``.
     * ``ypfin`` -- the forward solution derivative at the final time ``tB0``.
     * ``ySfin`` -- a pointer to an array of ``Ns`` vectors containing the sensitivities of the forward solution at the final time ``tB0``.
     * ``ypSfin`` -- a pointer to an array of ``Ns`` vectors containing the derivatives of the forward solution sensitivities at the final time ``tB0``.

   **Return value:**
     * ``IDA_NO_ADJ`` -- :c:func:`IDAAdjInit` has not been previously called.
     * ``IDA_ILL_INPUT`` -- Parameter ``which`` represented an invalid identifier, sensitivities were not active during forward integration, or :c:func:`IDAInitBS` or :c:func:`IDAReInitBS` has not been previously called.

   **Notes:**

   All failure return values are negative and therefore a test ``flag`` :math:`<
   0` will trap all :c:func:`IDACalcICBS` failures.  Note that
   :c:func:`IDACalcICBS` will correct the values of :math:`yB(tB_0)` and
   :math:`\dot{y}B(tB_0)` which were specified in the previous call to
   :c:func:`IDAInitBS` or :c:func:`IDAReInitBS`. To obtain the corrected values,
   call :c:func:`IDAGetConsistentICB` (see
   :numref:`IDAS.Usage.ADJ.user_callable.optional_ouput_b.iccalcB`).

   :c:func:`IDACalcICBS` will correct the values of :math:`yB(tB_0)` and
   :math:`\dot{y}B(tB_0)` which were specified in the previous call to
   :c:func:`IDAInitBS` or :c:func:`IDAReInitBS`. To obtain the corrected values,
   :call :c:func:`IDAGetConsistentICB`.


.. _IDAS.Usage.ADJ.user_callable.idasolveb:

Backward integration function
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The function :c:func:`IDASolveB` performs the integration of the backward problem. It
is essentially a wrapper for the IDAS main integration function :c:func:`IDASolve`
and, in the case in which checkpoints were needed, it evolves the solution of
the backward problem through a sequence of forward-backward integration pairs
between consecutive checkpoints. In each pair, the first run integrates the
original IVP forward in time and stores interpolation data; the second run
integrates the backward problem backward in time and performs the required
interpolation to provide the solution of the IVP to the backward problem.

The function :c:func:`IDASolveB` does not return the solution ``yB`` itself. To obtain
that, call the function :c:func:`IDAGetB`, which is also described below.

The :c:func:`IDASolveB` function does not support rootfinding, unlike :c:func:`IDASoveF`,
which supports the finding of roots of functions of :math:`(t,y,\dot{y})`. If
rootfinding was performed by :c:func:`IDASolveF`, then for the sake of efficiency, it
should be disabled for :c:func:`IDASolveB` by first calling :c:func:`IDARootInit` with
``nrtfn`` = 0.

The call to :c:func:`IDASolveB` has the form

.. c:function:: int IDASolveB(void * ida_mem, realtype tBout, int itaskB)

   The function :c:func:`IDASolveB` integrates the backward DAE problem.

   **Arguments:**
     * ``ida_mem`` -- pointer to the IDAS memory returned by :c:func:`IDACreate`.
     * ``tBout`` -- the next time at which a computed solution is desired.
     * ``itaskB`` --  output mode a flag indicating the job of the solver for the next step. The ``IDA_NORMAL`` task is to have the solver take internal steps until it has reached or just passed the user-specified value ``tBout``. The solver then interpolates in order to return an approximate value of :math:`yB(\texttt{tBout})`. The ``IDA_ONE_STEP`` option tells the solver to take just one internal step in the direction of ``tBout`` and return.

   **Return value:**
     * ``IDA_SUCCESS`` -- :c:func:`IDASolveB` succeeded.
     * ``IDA_MEM_NULL`` -- The ``ida_mem`` was ``NULL``.
     * ``IDA_NO_ADJ`` -- The function :c:func:`IDAAdjInit` has not been previously called.
     * ``IDA_NO_BCK`` -- No backward problem has been added to the list of backward problems by a call to :c:func:`IDACreateB`.
     * ``IDA_NO_FWD`` -- The function :c:func:`IDASolveF` has not been previously called.
     * ``IDA_ILL_INPUT`` -- One of the inputs to :c:func:`IDASolveB` is illegal.
     * ``IDA_BAD_ITASK`` -- The ``itaskB`` argument has an illegal value.
     * ``IDA_TOO_MUCH_WORK`` -- The solver took ``mxstep`` internal steps but could not reach ``tBout``.
     * ``IDA_TOO_MUCH_ACC`` -- The solver could not satisfy the accuracy demanded by the user for some internal step.
     * ``IDA_ERR_FAILURE`` -- Error test failures occurred too many times during one internal time step.
     * ``IDA_CONV_FAILURE`` -- Convergence test failures occurred too many times during one internal time step.
     * ``IDA_LSETUP_FAIL`` -- The linear solver's setup function failed in an unrecoverable manner.
     * ``IDA_SOLVE_FAIL`` -- The linear solver's solve function failed in an unrecoverable manner.
     * ``IDA_BCKMEM_NULL`` -- The IDAS memory for the backward problem was not created with a call to :c:func:`IDACreateB`.
     * ``IDA_BAD_TBOUT`` -- The desired output time ``tBout`` is outside the interval over which the forward problem was solved.
     * ``IDA_REIFWD_FAIL`` -- Reinitialization of the forward problem failed at the first checkpoint corresponding to the initial time of the forward problem.
     * ``IDA_FWD_FAIL`` -- An error occurred during the integration of the forward problem.

   **Notes:**
      All failure return values are negative and therefore a test
      ``flag``:math:`< 0`  will trap all :c:func:`IDASolveB` failures.  In the
      case of multiple checkpoints and multiple backward problems, a given  call
      to :c:func:`IDASolveB` in ``IDA_ONE_STEP`` mode may not advance every
      problem  one step, depending on the relative locations of the current
      times reached.  But repeated calls will eventually advance all problems to
      ``tBout``.


To obtain the solution ``yB`` to the backward problem, call the function
:c:func:`IDAGetB` as follows:

.. c:function:: int IDAGetB(void * ida_mem, int which, realtype * tret, N_Vector yB, N_Vector ypB)

   The function :c:func:`IDAGetB` provides the solution ``yB`` of the backward
   DAE  problem.

   **Arguments:**
     * ``ida_mem`` -- pointer to the IDAS memory returned by :c:func:`IDACreate`.
     * ``which`` -- the identifier of the backward problem.
     * ``tret`` -- the time reached by the solver output.
     * ``yB`` -- the backward solution at time ``tret``.
     * ``ypB`` -- the backward solution derivative at time ``tret``.

   **Return value:**
     * ``IDA_SUCCESS`` -- :c:func:`IDAGetB` was successful.
     * ``IDA_MEM_NULL`` -- ``ida_mem`` is ``NULL``.
     * ``IDA_NO_ADJ`` -- The function :c:func:`IDAAdjInit` has not been previously called.
     * ``IDA_ILL_INPUT`` -- The parameter ``which`` is an invalid identifier.

   **Notes:**
      To obtain the solution associated with a given backward problem at some
      other time within the last integration step, first obtain a pointer to the
      proper IDAS memory structure by calling :c:func:`IDAGetAdjIDABmem`  and
      then use it to call :c:func:`IDAGetDky`.

   .. warning::
      The user must allocate space for ``yB`` and ``ypB``.


.. _IDAS.Usage.ADJ.user_callable.optional_input_b:

Optional input functions for the backward problem
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

As for the forward problem there are numerous optional input parameters that
control the behavior of the IDAS solver for the backward problem. IDAS provides
functions that can be used to change these optional input parameters from their
default values which are then described in detail in the remainder of this
section, beginning with those for the main IDAS solver and continuing with those
for the linear solver interfaces. For the most casual use of IDAS, the reader
can skip to :numref:`IDAS.Usage.ADJ.user_supplied`.

We note that, on an error return, all of the optional input functions send an
error message to the error handler function. All error return values are
negative, so the test ``flag < 0`` will catch all errors. Finally, a call to a
``IDASet***B`` function can be made from the user’s calling program at any time
and, if successful, takes effect immediately.

Main solver optional input functions
""""""""""""""""""""""""""""""""""""

The adjoint module in IDAS provides wrappers for most of the optional input
functions defined in :numref:`IDAS.Usage.SIM.user_callable.optional_input`. The only
difference is that the user must specify the identifier ``which`` of the
backward problem within the list managed by IDAS.

The optional input functions defined for the backward problem are:

.. code-block:: c

     flag = IDASetUserDataB(ida_mem, which, user_dataB);
     flag = IDASetMaxOrdB(ida_mem, which, maxordB);
     flag = IDASetMaxNumStepsB(ida_mem, which, mxstepsB);
     flag = IDASetInitStepB(ida_mem, which, hinB)
     flag = IDASetMaxStepB(ida_mem, which, hmaxB);
     flag = IDASetSuppressAlgB(ida_mem, which, suppressalgB);
     flag = IDASetIdB(ida_mem, which, idB);
     flag = IDASetConstraintsB(ida_mem, which, constraintsB);

Their return value ``flag`` (of type ``int``) can have any of the return values
of their counterparts, but it can also be ``IDA_NO_ADJ`` if :c:func:`IDAAdjInit` has
not been called, or ``IDA_ILL_INPUT`` if ``which`` was an invalid identifier.

Linear solver interface optional input functions
""""""""""""""""""""""""""""""""""""""""""""""""

When using matrix-based linear solver modules for the backward problem, i.e., a
non-``NULL`` ``SUNMatrix`` object ``A`` was passed to :c:func:`IDASetLinearSolverB`,
the IDALS linear solver interface needs a function to compute an
approximation to the Jacobian matrix. This can be attached through a call to
either :c:func:`IDASetJacFnB` or :c:func:`IDASetJacFnBS`, with the second used when the
backward problem depends on the forward sensitivities.

.. c:function:: int IDASetJacFnB(void * ida_mem, int which, IDALsJacFnB jacB)

   The function :c:func:`IDASetJacFnB` specifies the Jacobian  approximation
   function to be used for the backward problem.

   **Arguments:**
     * ``ida_mem`` -- pointer to the IDAS memory block.
     * ``which`` -- represents the identifier of the backward problem.
     * ``jacB`` -- user-defined Jacobian approximation function.

   **Return value:**
     * ``IDALS_SUCCESS`` -- :c:func:`IDASetJacFnB` succeeded.
     * ``IDALS_MEM_NULL`` -- The ``ida_mem`` was ``NULL``.
     * ``IDALS_NO_ADJ`` -- The function :c:func:`IDAAdjInit` has not been previously called.
     * ``IDALS_LMEM_NULL`` -- The linear solver has not been initialized with a call to :c:func:`IDASetLinearSolverB`.
     * ``IDALS_ILL_INPUT`` -- The parameter ``which`` represented an invalid identifier.

   **Notes:**
      The previous routine ``IDADlsSetJacFnB`` is now a wrapper for this
      routine, and may still be used for backward-compatibility.  However,  this
      will be deprecated in future releases, so we recommend that  users
      transition to the new routine name soon.


.. c:function:: int IDASetJacFnBS(void * ida_mem, int which, IDALsJacFnBS jacBS)

   The function :c:func:`IDASetJacFnBS` specifies the Jacobian  approximation
   function to be used for the backward problem in the case  where the backward
   problem depends on the forward sensitivities.

   **Arguments:**
     * ``ida_mem`` -- pointer to the IDAS memory block.
     * ``which`` -- represents the identifier of the backward problem.
     * ``jacBS`` -- user-defined Jacobian approximation function.

   **Return value:**
     * ``IDALS_SUCCESS`` -- :c:func:`IDASetJacFnBS` succeeded.
     * ``IDALS_MEM_NULL`` -- The ``ida_mem`` was ``NULL``.
     * ``IDALS_NO_ADJ`` -- The function :c:func:`IDAAdjInit` has not been previously called.
     * ``IDALS_LMEM_NULL`` -- The linear solver has not been initialized with a call to :c:func:`IDASetLinearSolverBS`.
     * ``IDALS_ILL_INPUT`` -- The parameter ``which`` represented an invalid identifier.

   **Notes:**
      The previous routine, ``IDADlsSetJacFnBS``, is now deprecated.


The function :c:func:`IDASetLinearSolutionScalingB` can be used to enable or
disable solution scaling when using a matrix-based linear solver.

.. c:function:: int IDASetLinearSolutionScalingB(void * ida_mem, int which, booleantype onoffB)

   The function :c:func:`IDASetLinearSolutionScalingB` enables or disables
   scaling  the linear system solution to account for a change in :math:`\alpha`
   in the linear  system in the backward problem. For more details see :numref:`SUNLinSol.IDAS.Lagged`.

   **Arguments:**
     * ``ida_mem`` -- pointer to the IDAS memory block.
     * ``which`` -- represents the identifier of the backward problem.
     * ``onoffB`` -- flag to enable ``SUNTRUE`` or disable ``SUNFALSE`` scaling.

   **Return value:**
     * ``IDALS_SUCCESS`` -- The flag value has been successfully set.
     * ``IDALS_MEM_NULL`` -- The ``ida_mem`` pointer is ``NULL``.
     * ``IDALS_LMEM_NULL`` -- The IDALS linear solver interface has not been initialized.
     * ``IDALS_ILL_INPUT`` -- The attached linear solver is not matrix-based.

   **Notes:**

   By default scaling is enabled with matrix-based linear solvers when using
   BDF methods.

   By default scaling is enabled with matrix-based linear solvers when using BDF
   methods.

When using a matrix-free linear solver module for the backward problem,
the IDALS linear solver interface requires a function to compute an
approximation to the product between the Jacobian matrix :math:`J(t,y)` and a
vector :math:`v`. This may be performed internally using a difference-quotient
approximation, or it may be supplied by the user by calling one of the following
two functions:

.. c:function:: int IDASetJacTimesB(void * ida_mem, int which, \
                IDALsJacTimesSetupFnB jsetupB, IDALsJacTimesVecFnB jtimesB)

   The function :c:func:`IDASetJacTimesB` specifies the Jacobian-vector  setup
   and product functions to be used.

   **Arguments:**
     * ``ida_mem`` -- pointer to the IDAS memory block.
     * ``which`` -- the identifier of the backward problem.
     * ``jtsetupB`` -- user-defined function to set up the Jacobian-vector product. Pass ``NULL`` if no setup is necessary.
     * ``jtimesB`` -- user-defined Jacobian-vector product function.

   **Return value:**
     * ``IDALS_SUCCESS`` -- The optional value has been successfully set.
     * ``IDALS_MEM_NULL`` -- The ``ida_mem`` memory block pointer was ``NULL``.
     * ``IDALS_LMEM_NULL`` -- The IDALS linear solver has not been initialized.
     * ``IDALS_NO_ADJ`` -- The function :c:func:`IDAAdjInit` has not been previously called.
     * ``IDALS_ILL_INPUT`` -- The parameter ``which`` represented an invalid identifier.

   .. warning::

      The previous routine, ``IDASpilsSetJacTimesB``, is now deprecated.


.. c:function:: int IDASetJacTimesBS(void * ida_mem, int which, \
                IDALsJacTimesSetupFnBS jsetupBS, IDALsJacTimesVecFnBS jtimesBS)

   The function :c:func:`IDASetJacTimesBS` specifies the Jacobian-vector
   product setup and evaluation functions to be used, in the case where the
   backward problem depends on the forward sensitivities.

   **Arguments:**
     * ``ida_mem`` -- pointer to the IDAS memory block.
     * ``which`` -- the identifier of the backward problem.
     * ``jtsetupBS`` -- user-defined function to set up the Jacobian-vector product. Pass ``NULL`` if no setup is necessary.
     * ``jtimesBS`` -- user-defined Jacobian-vector product function.

   **Return value:**
     * ``IDALS_SUCCESS`` -- The optional value has been successfully set.
     * ``IDALS_MEM_NULL`` -- The ``ida_mem`` memory block pointer was ``NULL``.
     * ``IDALS_LMEM_NULL`` -- The IDALS linear solver has not been initialized.
     * ``IDALS_NO_ADJ`` -- The function :c:func:`IDAAdjInit` has not been previously called.
     * ``IDALS_ILL_INPUT`` -- The parameter ``which`` represented an invalid identifier.

   .. warning::

      The previous routine, ``IDASpilsSetJacTimesBS``, is now deprecated.


When using the default difference-quotient approximation to the Jacobian-vector
product for the backward problem, the user may specify the factor to use in
setting increments for the finite-difference approximation, via a call to
:c:func:`IDASetIncrementFactorB`.

.. c:function:: int IDASetIncrementFactorB(void * ida_mem, int which, realtype dqincfacB)

   The function :c:func:`IDASetIncrementFactorB` specifies the factor  in the
   increments used in the difference quotient approximations to matrix-vector
   products for the backward problem.  This routine can be used in both the
   cases where the backward problem  does and does not depend on the forward
   sensitvities.

   **Arguments:**
     * ``ida_mem`` -- pointer to the IDAS memory block.
     * ``which`` -- the identifier of the backward problem.
     * ``dqincfacB`` -- difference quotient approximation factor.

   **Return value:**
     * ``IDALS_SUCCESS`` -- The optional value has been successfully set.
     * ``IDALS_MEM_NULL`` -- The ``ida_mem`` pointer is ``NULL``.
     * ``IDALS_LMEM_NULL`` -- The IDALS linear solver has not been initialized.
     * ``IDALS_NO_ADJ`` -- The function :c:func:`IDAAdjInit` has not been previously called.
     * ``IDALS_ILL_INPUT`` -- The parameter ``which`` represented an invalid identifier.

   **Notes:**
      The default value is :math:`1.0`.

      The previous routine ``IDASpilsSetIncrementFactorB`` is now a deprecated.


Additionally, When using the internal difference quotient for the backward
problem, the user may also optionally supply an alternative residual function
for use in the Jacobian-vector product approximation by calling
:c:func:`IDASetJacTimesResFnB`. The alternative residual side function should
compute a suitable (and differentiable) approximation to the residual function
provided to :c:func:`IDAInitB` or :c:func:`IDAInitBS`. For example, as done in
:cite:p:`dorr2010numerical` for the forward integration of an ODE in explicit
form without sensitivity analysis, the alternative function may use lagged
values when evaluating a nonlinearity in the right-hand side to avoid
differencing a potentially non-differentiable factor.

.. c:function:: int IDASetJacTimesResFnB(void * ida_mem, int which, IDAResFn jtimesResFn)

   The function :c:func:`IDASetJacTimesResFnB` specifies an alternative DAE
   residual  function for use in the internal Jacobian-vector product difference
   quotient  approximation for the backward problem.

   **Arguments:**
     * ``ida_mem`` -- pointer to the IDAS memory block.
     * ``which`` -- the identifier of the backward problem.
     * ``jtimesResFn`` -- is the C function which computes the alternative DAE residual
       function to use in Jacobian-vector product difference quotient approximations. This
       function has the form ``res(t, yy, yp, resval, user_data)``. For full details see
       :numref:`IDAS.Usage.SIM.user_supplied.resFn`.

   **Return value:**
     * ``IDALS_SUCCESS`` -- The optional value has been successfully set.
     * ``IDALS_MEM_NULL`` -- The ``ida_mem`` pointer is ``NULL``.
     * ``IDALS_LMEM_NULL`` -- The IDALS linear solver has not been initialized.
     * ``IDALS_NO_ADJ`` -- The function :c:func:`IDAAdjInit` has not been previously called.
     * ``IDALS_ILL_INPUT`` -- The parameter ``which`` represented an invalid identifier or the internal difference quotient approximation is disabled.

   **Notes:**
      The default is to use the residual function provided to :c:func:`IDAInit`
      in the  internal difference quotient. If the input resudual function is
      ``NULL``,  the default is used.

      This function must be called *after* the
      IDALS linear solver interface  has been initialized through a call to
      :c:func:`IDASetLinearSolverB`.


When using an iterative linear solver for the backward problem, the user may
supply a preconditioning operator to aid in solution of the system, or she/he
may adjust the convergence tolerance factor for the iterative linear solver.
These may be accomplished through calling the following functions:

.. c:function:: int IDASetPreconditionerB(void * ida_mem, int which, IDALsPrecSetupFnB psetupB, IDALsPrecSolveFnB psolveB)

   The function :c:func:`IDASetPrecSolveFnB` specifies the preconditioner  setup
   and solve functions for the backward integration.

   **Arguments:**
     * ``ida_mem`` -- pointer to the IDAS memory block.
     * ``which`` -- the identifier of the backward problem.
     * ``psetupB`` -- user-defined preconditioner setup function.
     * ``psolveB`` -- user-defined preconditioner solve function.

   **Return value:**
     * ``IDALS_SUCCESS`` -- The optional value has been successfully set.
     * ``IDALS_MEM_NULL`` -- The ``ida_mem`` memory block pointer was ``NULL``.
     * ``IDALS_LMEM_NULL`` -- The IDALS linear solver has not been initialized.
     * ``IDALS_NO_ADJ`` -- The function :c:func:`IDAAdjInit` has not been previously called.
     * ``IDALS_ILL_INPUT`` -- The parameter ``which`` represented an invalid identifier.

   **Notes:**
      The ``psetupB`` argument may be ``NULL`` if no setup operation is involved
      in the preconditioner.

   .. warning::

      The previous routine ``IDASpilsSetPreconditionerB`` is now deprecated.


.. c:function:: int IDASetPreconditionerBS(void * ida_mem, int which, IDALsPrecSetupFnBS psetupBS, IDALsPrecSolveFnBS psolveBS)

   The function :c:func:`IDASetPrecSolveFnBS` specifies the preconditioner
   setup and solve functions for the backward integration, in the case  where
   the backward problem depends on the forward sensitivities.

   **Arguments:**
     * ``ida_mem`` -- pointer to the IDAS memory block.
     * ``which`` -- the identifier of the backward problem.
     * ``psetupBS`` -- user-defined preconditioner setup function.
     * ``psolveBS`` -- user-defined preconditioner solve function.

   **Return value:**
     * ``IDALS_SUCCESS`` -- The optional value has been successfully set.
     * ``IDALS_MEM_NULL`` -- The ``ida_mem`` memory block pointer was ``NULL``.
     * ``IDALS_LMEM_NULL`` -- The IDALS linear solver has not been initialized.
     * ``IDALS_NO_ADJ`` -- The function :c:func:`IDAAdjInit` has not been previously called.
     * ``IDALS_ILL_INPUT`` -- The parameter ``which`` represented an invalid identifier.

   **Notes:**
      The ``psetupBS`` argument may be ``NULL`` if no setup operation is
      involved  in the preconditioner.

   .. warning::

      The previous routine ``IDASpilsSetPreconditionerBS`` is now deprecated.


.. c:function:: int IDASetEpsLinB(void * ida_mem, int which, realtype eplifacB)

   The function :c:func:`IDASetEpsLinB` specifies the factor by  which the
   Krylov linear solver's convergence test constant is reduced  from the
   nonlinear iteration test constant. (See :numref:`IDAS.Mathematics.ivp_sol`).
   This routine can be used in both the cases wherethe backward problem does
   and does not depend on the forward sensitvities.

   **Arguments:**
     * ``ida_mem`` -- pointer to the IDAS memory block.
     * ``which`` -- the identifier of the backward problem.
     * ``eplifacB`` -- linear convergence safety factor :math:`>= 0.0`.

   **Return value:**
     * ``IDALS_SUCCESS`` -- The optional value has been successfully set.
     * ``IDALS_MEM_NULL`` -- The ``ida_mem`` pointer is ``NULL``.
     * ``IDALS_LMEM_NULL`` -- The IDALS linear solver has not been initialized.
     * ``IDALS_NO_ADJ`` -- The function :c:func:`IDAAdjInit` has not been previously called.
     * ``IDALS_ILL_INPUT`` -- The parameter ``which`` represented an invalid identifier.

   **Notes:**
      The default value is :math:`0.05`.

      Passing a value ``eplifacB`` :math:`= 0.0` also indicates using the
      default value.

   .. warning::

      The previous routine ``IDASpilsSetEpsLinB`` is now deprecated.

.. c:function:: int IDASetLSNormFactorB(void * ida_mem, int which, realtype nrmfac)

   The function :c:func:`IDASetLSNormFactorB` specifies the factor to use when
   converting from the integrator tolerance (WRMS norm) to the linear solver
   tolerance (L2 norm) for Newton linear system solves e.g.,  ``tol_L2 = fac *
   tol_WRMS``.  This routine can be used in both the cases wherethe backward
   problem  does and does not depend on the forward sensitvities.

   **Arguments:**
     * ``ida_mem`` -- pointer to the IDAS memory block.
     * ``which`` -- the identifier of the backward problem.
     * ``nrmfac`` -- the norm conversion factor. If ``nrmfac`` is:

       - :math:`> 0` then the provided value is used.
       - :math:`= 0` then the conversion factor is computed using the vector length i.e., ``nrmfac = N_VGetLength(y)`` default.
       - :math:`< 0` then the conversion factor is computed using the vector dot product ``nrmfac = N_VDotProd(v,v)`` where all the entries of ``v`` are one.

   **Return value:**
     * ``IDALS_SUCCESS`` -- The optional value has been successfully set.
     * ``IDALS_MEM_NULL`` -- The ``ida_mem`` pointer is ``NULL``.
     * ``IDALS_LMEM_NULL`` -- The IDALS linear solver has not been initialized.
     * ``IDALS_NO_ADJ`` -- The function :c:func:`IDAAdjInit` has not been previously called.
     * ``IDALS_ILL_INPUT`` -- The parameter ``which`` represented an invalid identifier.

   **Notes:**
      This function must be called after the IDALS linear solver  interface has
      been initialized through a call to  :c:func:`IDASetLinearSolverB`.

      Prior to the introduction of ``N_VGetLength`` in SUNDIALS v5.0.0 (IDAS
      v4.0.0) the value of ``nrmfac`` was computed using the vector dot product
      i.e., the ``nrmfac < 0`` case.


.. _IDAS.Usage.ADJ.user_callable.optional_output_b:

Optional output functions for the backward problem
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Main solver optional output functions
"""""""""""""""""""""""""""""""""""""

The user of the adjoint module in IDAS has access to any of the optional output
functions described in :numref:`IDAS.Usage.SIM.user_callable.optional_output`,
both for the main solver and for the linear solver modules. The first argument
of these ``IDAGet*`` and ``IDA*Get*`` functions is the pointer to the IDAS
memory block for the backward problem. In order to call any of these functions,
the user must first call the following function to obtain this pointer:

.. c:function:: int IDAGetAdjIDABmem(void * ida_mem, int which)

   The function :c:func:`IDAGetAdjIDABmem` returns a pointer to the IDAS  memory
   block for the backward problem.

   **Arguments:**
     * ``ida_mem`` -- pointer to the IDAS memory block created by :c:func:`IDACreate`.
     * ``which`` -- the identifier of the backward problem.

   **Return value:**
     * The return value, ``ida_memB`` (of type ``void *``), is a pointer to the
       idas memory for the backward problem.

   .. warning::
      The user should not modify ``ida_memB`` in any way.

      Optional output calls should pass ``ida_memB`` as the first argument;
      thus, for example, to get the number of integration steps: ``flag =
      IDAGetNumSteps(idas_memB,&nsteps)``.


To get values of the *forward* solution during a backward integration, use the
following function. The input value of ``t`` would typically be equal to that at
which the backward solution has just been obtained with :c:func:`IDAGetB`. In
any case, it must be within the last checkpoint interval used by
:c:func:`IDASolveB`.

.. c:function:: int IDAGetAdjY(void * ida_mem, realtype t, N_Vector y, N_Vector yp)

   The function :c:func:`IDAGetAdjY` returns the interpolated value of the
   forward solution :math:`y` and its derivative during a backward integration.

   **Arguments:**
     * ``ida_mem`` -- pointer to the IDAS memory block created by :c:func:`IDACreate`.
     * ``t`` -- value of the independent variable at which :math:`y` is desired input.
     * ``y`` -- forward solution :math:`y(t)`.
     * ``yp`` -- forward solution derivative :math:`\dot{y}(t)`.

   **Return value:**
     * ``IDA_SUCCESS`` -- :c:func:`IDAGetAdjY` was successful.
     * ``IDA_MEM_NULL`` -- ``ida_mem`` was ``NULL``.
     * ``IDA_GETY_BADT`` -- The value of ``t`` was outside the current checkpoint interval.

   .. warning::
      The user must allocate space for ``y`` and ``yp``.


.. c:function:: int IDAGetAdjCheckPointsInfo(void * ida_mem, IDAadjCheckPointRec *ckpnt)

   The function :c:func:`IDAGetAdjCheckPointsInfo` loads an array of
   ``ncheck+1``  records of type :c:func:`IDAadjCheckPointRec`.  The user must
   allocate space for the array ``ckpnt``.

   **Arguments:**
     * ``ida_mem`` -- pointer to the IDAS memory block created by :c:func:`IDACreate`.
     * ``ckpnt`` -- array of ``ncheck+1`` checkpoint records, each of type :c:func:`IDAadjCheckPointRec`.

   **Return value:**
     * ``void``

   **Notes:**
      The members of each record ``ckpnt[i]`` are:

      -  ``ckpnt[i].my_addr`` (``void *``) address of current checkpoint in ``ida_mem->ida_adj_mem``

      -  ``ckpnt[i].next_addr`` (``void *``) address of next checkpoint

      -  ``ckpnt[i].t0`` (``realtype``) start of checkpoint interval

      -  ``ckpnt[i].t1`` (``realtype``) end of checkpoint interval

      -  ``ckpnt[i].nstep`` (``long int``) step counter at ckeckpoint ``t0``

      -  ``ckpnt[i].order`` (``int``) method order at checkpoint ``t0``

      -  ``ckpnt[i].step`` (``realtype``) step size at checkpoint ``t0``


.. _IDAS.Usage.ADJ.user_callable.optional_ouput_b.iccalcB:

Initial condition calculation optional output function
"""""""""""""""""""""""""""""""""""""""""""""""""""""""

.. c:function:: int IDAGetConsistentICB(void * ida_mem, int which, N_Vector yB0_mod, N_Vector ypB0_mod)

   The function :c:func:`IDAGetConsistentICB` returns the corrected initial
   conditions  for backward problem calculated by :c:func:`IDACalcICB`.

   **Arguments:**
     * ``ida_mem`` -- pointer to the IDAS memory block.
     * ``which`` -- is the identifier of the backward problem.
     * ``yB0_mod`` -- consistent initial vector.
     * ``ypB0_mod`` -- consistent initial derivative vector.

   **Return value:**
     * IDA_SUCCESS -- The optional output value has been successfully set.
     * ``IDA_MEM_NULL`` -- The ``ida_mem`` pointer is ``NULL``.
     * ``IDA_NO_ADJ`` -- :c:func:`IDAAdjInit` has not been previously called.
     * ``IDA_ILL_INPUT`` -- Parameter ``which`` did not refer a valid backward problem identifier.

   **Notes:**
      If the consistent solution vector or consistent derivative vector  is not
      desired, pass ``NULL`` for the corresponding argument.

      .. warning::
         The user must allocate space for ``yB0_mod`` and ``ypB0_mod``  (if not
         ``NULL``).


.. _IDAS.Usage.ADJ.user_callable.backquad:

Backward integration of quadrature equations
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Not only the backward problem but also the backward quadrature equations may or
may not depend on the forward sensitivities. Accordingly, one of the
:c:func:`IDAQuadInitB` or :c:func:`IDAQuadInitBS` should be used to allocate
internal memory and to initialize backward quadratures. For any other operation
(extraction, optional input/output, reinitialization, deallocation), the same
function is called regardless of whether or not the quadratures are
sensitivity-dependent.

.. _IDAS.Usage.ADJ.user_callable.backquad.idaquadinitb:

Backward quadrature initialization functions
""""""""""""""""""""""""""""""""""""""""""""

The function :c:func:`IDAQuadInitB` initializes and allocates memory for the
backward integration of quadrature equations that do not depende on forward
sensititvities. It has the following form:

.. c:function:: int IDAQuadInitB(void * ida_mem, int which, IDAQuadRhsFnB rhsQB, N_Vector yQB0)

   The function :c:func:`IDAQuadInitB` provides required problem specifications,
   allocates internal memory, and initializes backward quadrature integration.

   **Arguments:**
     * ``ida_mem`` -- pointer to the IDAS memory block.
     * ``which`` -- the identifier of the backward problem.
     * ``rhsQB`` -- is the C function which computes :math:`fQB` , the residual of the backward quadrature
       equations. This function has the form ``rhsQB(t, y, yp, yB, ypB, rhsvalBQ, user_dataB)`` see
       :numref:`IDAS.Usage.ADJ.RHS_quad_B`.
     * ``yQB0`` -- is the value of the quadrature variables at ``tB0``.

   **Return value:**
     * ``IDA_SUCCESS`` -- The call to :c:func:`IDAQuadInitB` was successful.
     * ``IDA_MEM_NULL`` -- ``ida_mem`` was ``NULL``.
     * ``IDA_NO_ADJ`` -- The function :c:func:`IDAAdjInit` has not been previously called.
     * ``IDA_MEM_FAIL`` -- A memory allocation request has failed.
     * ``IDA_ILL_INPUT`` -- The parameter ``which`` is an invalid identifier.

.. c:function:: int IDAQuadInitBS(void * ida_mem, int which, IDAQuadRhsFnBS rhsQBS, N_Vector yQBS0)

   The function :c:func:`IDAQuadInitBS` provides required problem
   specifications,  allocates internal memory, and initializes backward
   quadrature integration with sensitivities.

   **Arguments:**
     * ``ida_mem`` -- pointer to the IDAS memory block.
     * ``which`` -- the identifier of the backward problem.
     * ``rhsQBS`` -- is the C function which computes :math:`fQBS`, the residual of the backward quadrature
       equations. This function has the form ``rhsQBS(t, y, yp, yS, ypS, yB, ypB, rhsvalBQS, user_dataB)``
       see :numref:`IDAS.Usage.ADJ.RHS_quad_sens_B`.
     * ``yQBS0`` -- is the value of the sensitivity-dependent quadrature variables at ``tB0``.

   **Return value:**
     * ``IDA_SUCCESS`` -- The call to :c:func:`IDAQuadInitBS` was successful.
     * ``IDA_MEM_NULL`` -- ``ida_mem`` was ``NULL``.
     * ``IDA_NO_ADJ`` -- The function :c:func:`IDAAdjInit` has not been previously called.
     * ``IDA_MEM_FAIL`` -- A memory allocation request has failed.
     * ``IDA_ILL_INPUT`` -- The parameter ``which`` is an invalid identifier.


The integration of quadrature equations during the backward phase can be
re-initialized by calling the following function. Before calling
:c:func:`IDAQuadReInitB` for a new backward problem, call any desired solution
extraction functions ``IDAGet**`` associated with the previous backward problem.

.. c:function:: int IDAQuadReInitB(void * ida_mem, int which, N_Vector yQB0)

   The function :c:func:`IDAQuadReInitB` re-initializes the backward quadrature integration.

   **Arguments:**
     * ``ida_mem`` -- pointer to the IDAS memory block.
     * ``which`` -- the identifier of the backward problem.
     * ``yQB0`` -- is the value of the quadrature variables at ``tB0``.

   **Return value:**
     * ``IDA_SUCCESS`` -- The call to :c:func:`IDAQuadReInitB` was successful.
     * ``IDA_MEM_NULL`` -- ``ida_mem`` was ``NULL``.
     * ``IDA_NO_ADJ`` -- The function :c:func:`IDAAdjInit` has not been previously called.
     * ``IDA_MEM_FAIL`` -- A memory allocation request has failed.
     * ``IDA_NO_QUAD`` -- Quadrature integration was not activated through a previous call to :c:func:`IDAQuadInitB`.
     * ``IDA_ILL_INPUT`` -- The parameter ``which`` is an invalid identifier.

   **Notes:**
      :c:func:`IDAQuadReInitB` can be used after a call to either :c:func:`IDAQuadInitB`  or :c:func:`IDAQuadInitBS`.


.. _IDAS.Usage.ADJ.user_callable.backquad.quad_get_b:

Backward quadrature extraction function
"""""""""""""""""""""""""""""""""""""""

To extract the values of the quadrature variables at the last return time of
:c:func:`IDASolveB`, IDAS provides a wrapper for the function
:c:func:`IDAGetQuad`. The call to this function has the form


.. c:function:: int IDAGetQuadB(void * ida_mem, int which, realtype * tret, N_Vector yQB)

   The function :c:func:`IDAGetQuadB` returns the quadrature solution vector
   after  a successful return from :c:func:`IDASolveB`.

   **Arguments:**
     * ``ida_mem`` -- pointer to the IDAS memory.
     * ``tret`` -- the time reached by the solver output.
     * ``which`` -- the identifier of the backward problem.
     * ``yQB`` -- the computed quadrature vector.

   **Return value:**
     * ``IDA_SUCCESS`` -- :c:func:`IDAGetQuadB` was successful.
     * ``IDA_MEM_NULL`` -- ``ida_mem`` is ``NULL``.
     * ``IDA_NO_ADJ`` -- The function :c:func:`IDAAdjInit` has not been previously called.
     * IDA_NO_QUAD -- Quadrature integration was not initialized.
     * IDA_BAD_DKY -- ``yQB`` was ``NULL``.
     * ``IDA_ILL_INPUT`` -- The parameter ``which`` is an invalid identifier.

   **Notes:**
      To obtain the quadratures associated with a given backward problem at some
      other time within the last integration step, first obtain a pointer to the
      proper IDAS memory structure by calling :c:func:`IDAGetAdjIDABmem`  and
      then use it to call :c:func:`IDAGetQuadDky`.

      .. warning::
        The user must allocate space for ``yQB``.


.. _IDAS.Usage.ADJ.user_callable.backquad.quad_optional_input_B:

Optional input/output functions for backward quadrature integration
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Optional values controlling the backward integration of quadrature equations can
be changed from their default values through calls to one of the following
functions which are wrappers for the corresponding optional input functions
defined in :numref:`IDAS.Usage.Purequad.quad_optional_input`. The user
must specify the identifier ``which`` of the backward problem for which the
optional values are specified.

.. code-block:: c

     flag = IDASetQuadErrConB(ida_mem, which, errconQ);
     flag = IDAQuadSStolerancesB(ida_mem, which, reltolQ, abstolQ);
     flag = IDAQuadSVtolerancesB(ida_mem, which, reltolQ, abstolQ);

Their return value ``flag`` (of type ``int``) can have any of the return values
of its counterparts, but it can also be ``IDA_NO_ADJ`` if the function
:c:func:`IDAAdjInit` has not been previously called or ``IDA_ILL_INPUT`` if the
parameter ``which`` was an invalid identifier.

Access to optional outputs related to backward quadrature integration can be
obtained by calling the corresponding ``IDAGetQuad*`` functions (see
:numref:`IDAS.Usage.Purequad.quad_optional_output`). A pointer ``ida_memB`` to
the IDAS memory block for the backward problem, required as the first argument
of these functions, can be obtained through a call to the functions
:c:func:`IDAGetAdjIDABmem`.


.. _IDAS.Usage.ADJ.user_supplied:

User-supplied functions for adjoint sensitivity analysis
--------------------------------------------------------

In addition to the required DAE residual function and any optional functions for
the forward problem, when using the adjoint sensitivity module in IDAS, the user
must supply one function defining the backward problem DAE and, optionally,
functions to supply Jacobian-related information and one or two functions that
define the preconditioner (if applicable for the choice of ``SUNLinearSolver``
object) for the backward problem. Type definitions for all these user-supplied
functions are given below.

.. _IDAS.Usage.ADJ.user_supplied.DAEres_b:

DAE residual for the backward problem
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The user must provide a ``resB`` function of type ``IDAResFnB`` defined as follows:

.. c:type:: int (*IDAResFnB)(realtype t, N_Vector y, N_Vector yp, N_Vector yB, N_Vector ypB, N_Vector resvalB, void *user_dataB)

   This function evaluates the residual of the backward problem DAE system.
   This could be :eq:`IDAS_adj_eqns` or :eq:`IDAS_adj1_eqns`.

   **Arguments:**
     * ``t`` -- is the current value of the independent variable.
     * ``y`` -- is the current value of the forward solution vector.
     * ``yp`` -- is the current value of the forward solution derivative vector.
     * ``yB`` -- is the current value of the backward dependent variable vector.
     * ``ypB`` -- is the current value of the backward dependent derivative vector.
     * ``resvalB`` -- is the output vector containing the residual for the backward DAE problem.
     * ``user_dataB`` -- is a pointer to user data, same as passed to :c:func:`IDASetUserDataB` .

   **Return value:**
      An ``IDAResFnB`` should return 0 if successful, a positive value if a recoverable
      error occurred (in which case IDAS will attempt to correct), or a negative
      value if an unrecoverabl failure occurred (in which case the integration stops and
      :c:func:`IDASolveB` returns ``IDA_RESFUNC_FAIL``).

   **Notes:**
      Allocation of memory for ``resvalB`` is handled within IDAS.  The ``y``,
      ``yp``, ``yB``, ``ypB``, and ``resvalB`` arguments are all  of type
      ``N_Vector``, but ``yB``, ``ypB``, and ``resvalB`` typically have
      different internal representations from ``y`` and ``yp``.  It is the
      user's  responsibility to access the vector data consistently (including
      the use of the  correct accessor macros from each ``N_Vector``
      implementation). The ``user_dataB`` pointer is passed to the user's ``resB``
      function every time it is called and can be the same as the  ``user_data``
      pointer used for the forward problem.

      .. warning::
        Before calling the user's ``resB`` function, IDAS needs to evaluate
        (through interpolation) the values of the states from the forward
        integration.  If an error occurs in the interpolation, IDAS triggers an
        unrecoverable  failure in the residual function which will halt the
        integration and  :c:func:`IDASolveB` will return ``IDA_RESFUNC_FAIL``.


.. _IDAS.Usage.ADJ.DAEres_bs:

DAE residual for the backward problem depending on the forward sensitivities
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The user must provide a ``resBS`` function of type ``IDAResFnBS`` defined as
follows:

.. c:type:: int (*IDAResFnBS)(realtype t, N_Vector y, N_Vector yp, N_Vector *yS, N_Vector *ypS, N_Vector yB, N_Vector ypB, N_Vector resvalB, void *user_dataB)

   This function evaluates the residual of the backward problem DAE system.
   This could be :eq:`IDAS_adj_eqns` or :eq:`IDAS_adj1_eqns`.

   **Arguments:**
     * ``t`` -- is the current value of the independent variable.
     * ``y`` -- is the current value of the forward solution vector.
     * ``yp`` -- is the current value of the forward solution derivative vector.
     * ``yS`` -- a pointer to an array of ``Ns`` vectors containing the sensitivities of    the forward solution.
     * ``ypS`` -- a pointer to an array of ``Ns`` vectors containing the derivatives of    the forward sensitivities.
     * ``yB`` -- is the current value of the backward dependent variable vector.
     * ``ypB`` -- is the current value of the backward dependent derivative vector.
     * ``resvalB`` -- is the output vector containing the residual for the backward DAE problem.
     * ``user_dataB`` -- is a pointer to user data, same as passed to :c:func:`IDASetUserDataB` .

   **Return value:**
      An ``IDAResFnBS`` should return 0 if successful, a positive value if a
      recoverable error occurred (in which case IDAS will attempt to correct),
      or a negative value if an unrecoverable error occurred (in which case the
      integration stops and :c:func:`IDASolveB` returns ``IDA_RESFUNC_FAIL``).

   **Notes:**
      Allocation of memory for ``resvalB`` is handled within IDAS.  The ``y``,
      ``yp``, ``yB``, ``ypB``, and ``resvalB`` arguments are all  of type
      ``N_Vector``, but ``yB``, ``ypB``, and ``resvalB`` typically have
      different internal representations from ``y`` and ``yp``. Likewise for
      each  ``yS[i]`` and ``ypS[i]``. It is the user's  responsibility to access
      the vector data consistently (including the use of the  correct accessor
      macros from each ``N_Vector`` implementation).  The ``user_dataB`` pointer
      is passed to  the user's ``resBS`` function every time it is called and
      can be the same as the  ``user_data`` pointer used for the forward
      problem.

      .. warning::
        Before calling the user's ``resBS`` function, IDAS needs to evaluate
        (through interpolation) the values of the states from the forward
        integration.  If an error occurs in the interpolation, IDAS triggers an
        unrecoverable  failure in the residual function which will halt the
        integration and  :c:func:`IDASolveB` will return ``IDA_RESFUNC_FAIL``.


.. _IDAS.Usage.ADJ.RHS_quad_B:

Quadrature right-hand side for the backward problem
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The user must provide an ``fQB`` function of type ``IDAQuadRhsFnB`` defined by

.. c:type:: int (*IDAQuadRhsFnB)(realtype t, N_Vector y, N_Vector yp, N_Vector yB, N_Vector ypB, N_Vector rhsvalBQ, void *user_dataB)

   This function computes the quadrature equation right-hand side for the
   backward problem.

   **Arguments:**
     * ``t`` -- is the current value of the independent variable.
     * ``y`` -- is the current value of the forward solution vector.
     * ``yp`` -- is the current value of the forward solution derivative vector.
     * ``yB`` -- is the current value of the backward dependent variable vector.
     * ``ypB`` -- is the current value of the backward dependent derivative vector.
     * ``rhsvalBQ`` -- is the output vector containing the residual for the backward quadrature    equations.
     * ``user_dataB`` -- is a pointer to user data, same as passed to :c:func:`IDASetUserDataB` .

   **Return value:**
      An ``IDAQuadRhsFnB`` should return 0 if successful, a positive value if a
      recoverable error occurred (in which case IDAS will attempt to correct),
      or a negative value if it failed unrecoverably (in which case the
      integration is halted and :c:func:`IDASolveB` returns
      ``IDA_QRHSFUNC_FAIL``).

   **Notes:**
      Allocation of memory for ``rhsvalBQ`` is handled within IDAS.  The ``y``,
      ``yp``, ``yB``, ``ypB``, and ``rhsvalBQ`` arguments are all  of type
      ``N_Vector``, but they typically all have  different internal
      representations. It is the user's  responsibility to access the vector
      data consistently (including the use of the  correct accessor macros from
      each ``N_Vector`` implementation). For the sake of  computational
      efficiency, the vector functions in the two ``N_Vector`` implementations
      provided with IDAS do not perform any consistency checks with repsect to
      their  ``N_Vector`` arguments (see :numref:`NVectors`).  The ``user_dataB``
      pointer is passed to the user's ``fQB`` function every time  it is called
      and can be the same as the ``user_data`` pointer used for the forward
      problem.

      .. warning::
        Before calling the user's ``fQB`` function, IDAS needs to evaluate
        (through interpolation) the values of the states from the forward
        integration.  If an error occurs in the interpolation, IDAS triggers an
        unrecoverable  failure in the quadrature right-hand side function which
        will halt the integration and  :c:func:`IDASolveB` will return
        ``IDA_QRHSFUNC_FAIL``.


.. _IDAS.Usage.ADJ.RHS_quad_sens_B:

Sensitivity-dependent quadrature right-hand side for the backward problem
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The user must provide an ``fQBS`` function of type ``IDAQuadRhsFnBS`` defined by

.. c:type:: int (*IDAQuadRhsFnBS)(realtype t, N_Vector y, N_Vector yp, N_Vector *yS, N_Vector *ypS, N_Vector yB, N_Vector ypB, N_Vector rhsvalBQS, void *user_dataB)

   This function computes the quadrature equation residual for the  backward problem.

   **Arguments:**
     * ``t`` -- is the current value of the independent variable.
     * ``y`` -- is the current value of the forward solution vector.
     * ``yp`` -- is the current value of the forward solution derivative vector.
     * ``yS`` -- a pointer to an array of ``Ns`` vectors containing the sensitivities of    the forward solution.
     * ``ypS`` -- a pointer to an array of ``Ns`` vectors containing the derivatives of    the forward sensitivities.
     * ``yB`` -- is the current value of the backward dependent variable vector.
     * ``ypB`` -- is the current value of the backward dependent derivative vector.
     * ``rhsvalBQS`` -- is the output vector containing the residual for the backward quadrature    equations.
     * ``user_dataB`` -- is a pointer to user data, same as passed to :c:func:`IDASetUserDataB` .

   **Return value:**
      An ``IDAQuadRhsFnBS`` should return 0 if successful, a positive value if a
      recoverable error occurred (in which case IDAS will attempt to correct),
      or a negative value if it failed unrecoverably (in which case the
      integration is halted and :c:func:`IDASolveB` returns
      ``IDA_QRHSFUNC_FAIL``).

   **Notes:**
      Allocation of memory for ``rhsvalBQS`` is handled within IDAS.  The ``y``,
      ``yp``, ``yB``, ``ypB``, and ``rhsvalBQS`` arguments are all  of type
      ``N_Vector``, but they typically do not all have the same internal
      representations.  Likewise for each ``yS[i]`` and ``ypS[i]``.  It is the
      user's  responsibility to access the vector data consistently (including
      the use of the  correct accessor macros from each ``N_Vector``
      implementation).  The ``user_dataB`` pointer is passed to
      the user's ``fQBS`` function every time  it is called and can be the same
      as the ``user_data`` pointer used for the forward  problem.

      .. warning::
        Before calling the user's ``fQBS`` function, IDAS needs to evaluate
        (through interpolation) the values of the states from the forward
        integration.  If an error occurs in the interpolation, IDAS triggers an
        unrecoverable  failure in the quadrature right-hand side function which
        will halt the integration and  :c:func:`IDASolveB` will return
        ``IDA_QRHSFUNC_FAIL``.


.. _IDAS.Usage.ADJ.jacFn_b:

Jacobian construction for the backward problem (matrix-based linear solvers)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If a matrix-based linear solver module is is used for the backward problem
(i.e., :c:func:`IDASetLinearSolverB` is called with non-``NULL`` ``SUNMatrix``
argument in the step described in :numref:`IDAS.Usage.ADJ.skeleton_adj`), the
user may provide a function of type ``IDALsJacFnB`` or :c:type:`IDALsJacFnBS`, defined
as follows:

.. c:type:: int (*IDALsJacFnB)(realtype tt, realtype c_jB, N_Vector yy, N_Vector yp, N_Vector yyB, N_Vector ypB, N_Vector rrB, SUNMatrix JacB, void *user_dataB, N_Vector tmp1B, N_Vector tmp2B, N_Vector tmp3B)

   This function computes the Jacobian of the backward problem (or an
   approximation  to it).

   **Arguments:**
     * ``tt`` -- is the current value of the independent variable.
     * ``c_jB`` -- is the scalar in the system Jacobian, proportional to the inverse of the  step size (:math:`\alpha` in :eq:`IDAS_DAE_Jacobian`).
     * ``yy`` -- is the current value of the forward solution vector.
     * ``yp`` -- is the current value of the forward solution derivative vector.
     * ``yB`` -- is the current value of the backward dependent variable vector.
     * ``ypB`` -- is the current value of the backward dependent derivative vector.
     * ``rrB`` -- is the current value of the residual for the backward problem.
     * ``JacB`` -- is the output approximate Jacobian matrix.
     * ``user_dataB`` -- is a pointer to user data — the parameter passed to :c:func:`IDASetUserDataB` .
     * ``tmp1B``, ``tmp2B``, ``tmp3B`` -- are pointers to memory allocated for variables of type ``N_Vector`` which can be used by the :c:type:`IDALsJacFnB` function    as temporary storage or work space.

   **Return value:**
      An :c:type:`IDALsJacFnB` should return 0 if successful, a positive value if a
      recoverable error occurred (in which case IDAS will attempt to correct,
      while IDALS sets ``last_flag`` to ``IDALS_JACFUNC_RECVR``), or a negative
      value if it failed unrecoverably (in which case the integration is halted,
      :c:func:`IDASolveB` returns ``IDA_LSETUP_FAIL`` and IDALS sets
      ``last_flag`` to ``IDALS_JACFUNC_UNRECVR``).

   **Notes:**
      A user-supplied Jacobian function must load the  matrix ``JacB`` with an
      approximation to the Jacobian matrix  at the point ``(tt, yy, yB)``,
      where ``yy`` is the solution  of the original IVP at time ``tt``, and
      ``yB`` is the solution of the  backward problem at the same time.
      Information regarding the structure of the specific ``SUNMatrix``
      structure (e.g. number of rows, upper/lower bandwidth, sparsity  type) may
      be obtained through using the implementation-specific  ``SUNMatrix``
      interface functions (see Chapter :numref:`SUNMatrix` for  details).  With direct linear
      solvers (i.e., linear solvers with type  ``SUNLINEARSOLVER_DIRECT``), the
      Jacobian matrix :math:`J(t,y)` is zeroed out  prior to calling the
      user-supplied Jacobian function so only nonzero elements  need to be
      loaded into ``JacB``.

      .. warning::
        Before calling the user's ``IDALsJacFnB``, IDAS needs to evaluate
        (through interpolation) the values of the states from the forward
        integration.  If an error occurs in the interpolation, IDAS triggers an
        unrecoverable  failure in the Jacobian function which will halt the
        integration  (:c:func:`IDASolveB` returns ``IDA_LSETUP_FAIL`` and IDALS
        sets ``last_flag`` to  ``IDALS_JACFUNC_UNRECVR``).

        The previous
        function type ``IDADlsJacFnB`` is identical to  ``IDALsJacFnB``, and may
        still be used for backward-compatibility.  However, this will be
        deprecated in future releases, so we recommend  that users transition to
        the new function type name soon.


.. c:type:: int (*IDALsJacFnBS)(realtype tt, realtype c_jB, N_Vector yy, N_Vector yp, N_Vector *yS, N_Vector *ypS, N_Vector yyB, N_Vector ypB, N_Vector rrB, SUNMatrix JacB, void *user_dataB, N_Vector tmp1B, N_Vector tmp2B, N_Vector tmp3B);

   This function computes the Jacobian of the backward problem (or an
   approximation to it), in the case where the backward problem depends on the
   forward sensitivities.

   **Arguments:**
     * ``tt`` -- is the current value of the independent variable.
     * ``c_jB`` -- is the scalar in the system Jacobian, proportional to the inverse of the step size (:math:`\alpha` in :eq:`IDAS_DAE_Jacobian`).
     * ``yy`` -- is the current value of the forward solution vector.
     * ``yp`` -- is the current value of the forward solution derivative vector.
     * ``yS`` -- a pointer to an array of ``Ns`` vectors containing the sensitivities    of the forward solution.
     * ``ypS`` -- a pointer to an array of ``Ns`` vectors containing the derivatives    of the forward solution sensitivities.
     * ``yB`` -- is the current value of the backward dependent variable vector.
     * ``ypB`` -- is the current value of the backward dependent derivative vector.
     * ``rrb`` -- is the current value of the residual for the backward problem.
     * ``JacB`` -- is the output approximate Jacobian matrix.
     * ``user_dataB`` -- is a pointer to user data — the parameter passed to :c:func:`IDASetUserDataB` .
     * ``tmp1B``, ``tmp2B``, ``tmp3B`` -- are pointers to memory allocated  for variables of type ``N_Vector`` which    can be used by :c:type:`IDALsJacFnBS` as temporary storage or work space.

   **Return value:**
      An :c:type:`IDALsJacFnBS` should return 0 if successful, a positive value if a
      recoverable error occurred (in which case IDAS will attempt to correct,
      while IDALS sets ``last_flag`` to ``IDALS_JACFUNC_RECVR``), or a negative
      value if it failed unrecoverably (in which case the integration is halted,
      :c:func:`IDASolveB` returns ``IDA_LSETUP_FAIL`` and IDALS sets
      ``last_flag`` to ``IDALS_JACFUNC_UNRECVR``).

   **Notes:**
      A user-supplied dense Jacobian function must load the  matrix ``JacB``
      with an approximation to the Jacobian matrix  at the point
      ``(tt, yy, yS, yB)``, where ``yy`` is the solution  of the
      original IVP at time ``tt``, ``yS`` is the array of forward sensitivities
      at time ``tt``, and ``yB`` is the solution of the backward problem at the
      same time.  Information regarding the structure of the specific
      ``SUNMatrix``  structure (e.g. number of rows, upper/lower bandwidth,
      sparsity  type) may be obtained through using the implementation-specific
      ``SUNMatrix`` interface functions (see Chapter :numref:`SUNMatrix` for  details).  With
      direct linear solvers (i.e., linear solvers with type
      ``SUNLINEARSOLVER_DIRECT``, the Jacobian matrix :math:`J(t,y)` is zeroed
      out prior  to calling the user-supplied Jacobian function so only nonzero
      elements need  to be loaded into ``JacB``.

      .. warning::
        Before calling the user's :c:type:`IDALsJacFnBS`, IDAS needs to evaluate
        (through interpolation) the values of the states from the forward
        integration.  If an error occurs in the interpolation, IDAS triggers an
        unrecoverable  failure in the Jacobian function which will halt the
        integration  (:c:func:`IDASolveB` returns ``IDA_LSETUP_FAIL`` and IDALS
        sets ``last_flag`` to  ``IDALS_JACFUNC_UNRECVR``).

        The previous
        function type ``IDADlsJacFnBS`` is identical to  :c:type:`IDALsJacFnBS`, and
        may still be used for backward-compatibility.  However, this will be
        deprecated in future releases, so we recommend  that users transition to
        the new function type name soon.


.. _IDAS.Usage.ADJ.user_supplied.jactimesvec_b:

Jacobian-vector product for the backward problem (matrix-free linear solvers)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If a matrix-free linear solver is selected for the backward problem (i.e.,
:c:func:`IDASetLinearSolverB` is called with ``NULL``-valued ``SUNMatrix``
argument in the steps described in :numref:`IDAS.Usage.ADJ.skeleton_adj`), the user may
provide a function of type ``IDALsJacTimesVecFnB`` or ``IDALsJacTimesVecFnBS``
in the following form, to compute matrix-vector products :math:`Jv`. If such a
function is not supplied, the default is a difference quotient approximation to
these products.

.. c:type:: int (*IDALsJacTimesVecFnB)(realtype t, N_Vector yy, N_Vector yp, N_Vector yB, N_Vector ypB, N_Vector resvalB, N_Vector vB, N_Vector JvB, realtype cjB, void *user_dataB, N_Vector tmp1B, N_Vector tmp2B)

   This function computes the action of the backward problem Jacobian ``JB``  on
   a given vector ``vB``.

   **Arguments:**
     * ``t`` -- is the current value of the independent variable.
     * ``yy`` -- is the current value of the forward solution vector.
     * ``yp`` -- is the current value of the forward solution derivative vector.
     * ``yB`` -- is the current value of the backward dependent variable vector.
     * ``ypB`` -- is the current value of the backward dependent derivative vector.
     * ``resvalB`` -- is the current value of the residual for the backward problem.
     * ``vB`` -- is the vector by which the Jacobian must be multiplied.
     * ``JvB`` -- is the computed output vector, ``JB*vB`` .
     * ``cjB`` -- is the scalar in the system Jacobian, proportional to the inverse of the    step size ( :math:`\alpha` in :eq:`IDAS_DAE_Jacobian` ).
     * ``user_dataB`` -- is a pointer to user data — the same as the ``user_dataB`` parameter passed to :c:func:`IDASetUserDataB` .
     * ``tmp1B``, ``tmp2B`` -- are pointers to memory allocated for variables of type ``N_Vector`` which    can be used by ``IDALsJacTimesVecFnB`` as temporary storage or work space.

   **Return value:**
      The return value of a function of type ``IDALsJtimesVecFnB`` should be if
      successful or nonzero if an error was encountered, in which case the
      integration is halted.

   **Notes:**
      A user-supplied Jacobian-vector product function must load the vector
      ``JvB``  with the product of the Jacobian of the backward problem  at the
      point ``(t, y, yB)`` and the vector ``vB``.  Here, ``y`` is the
      solution of the original IVP at time ``t`` and  ``yB`` is the solution of
      the backward problem at the same time.  The rest of the arguments are
      equivalent to those passed to a function of type  ``IDALsJacTimesVecFn``
      (see :numref:`IDAS.Usage.SIM.user_supplied.jtimesFn`).  If the backward
      problem is the adjoint of :math:`{\dot y} = f(t, y)`, then this function
      is to compute :math:`-\left({\partial f}/{\partial y_i}\right)^T v_B`.

   .. warning::

      The previous function type ``IDASpilsJacTimesVecFnB`` is identical to
      ``IDALsJacTimesVecFnB``, and may still be used for
      backward-compatibility.  However, this will be deprecated in future
      releases, so we recommend that users transition to the new function  type
      name soon.


.. c:type:: int (*IDALsJacTimesVecFnBS)(realtype t, N_Vector yy, N_Vector yp, N_Vector *yyS, N_Vector *ypS, N_Vector yB, N_Vector ypB, N_Vector resvalB, N_Vector vB, N_Vector JvB, realtype cjB, void *user_dataB, N_Vector tmp1B, N_Vector tmp2B)

   This function computes the action of the backward problem Jacobian ``JB``  on
   a given vector ``vB``, in the case where the backward problem depends  on the
   forward sensitivities.

   **Arguments:**
     * ``t`` -- is the current value of the independent variable.
     * ``yy`` -- is the current value of the forward solution vector.
     * ``yp`` -- is the current value of the forward solution derivative vector.
     * ``yyS`` -- a pointer to an array of ``Ns`` vectors containing the sensitivities of    the forward solution.
     * ``ypS`` -- a pointer to an array of ``Ns`` vectors containing the derivatives of    the forward sensitivities.
     * ``yB`` -- is the current value of the backward dependent variable vector.
     * ``ypB`` -- is the current value of the backward dependent derivative vector.
     * ``resvalB`` -- is the current value of the residual for the backward problem.
     * ``vB`` -- is the vector by which the Jacobian must be multiplied.
     * ``JvB`` -- is the computed output vector, ``JB*vB`` .
     * ``cjB`` -- is the scalar in the system Jacobian, proportional to the inverse of the    step size ( :math:`\alpha` in :eq:`IDAS_DAE_Jacobian` ).
     * ``user_dataB`` -- is a pointer to user data — the same as the ``user_dataB`` parameter passed to :c:func:`IDASetUserDataB` .
     * ``tmp1B``, ``tmp2B`` -- are pointers to memory allocated for variables of type ``N_Vector`` which    can be used by ``IDALsJacTimesVecFnBS`` as temporary storage or work space.

   **Return value:**
      The return value of a function of type ``IDALsJtimesVecFnBS`` should be if
      successful or nonzero if an error was encountered, in which case the
      integration is halted.

   **Notes:**
      A user-supplied Jacobian-vector product function must load the vector
      ``JvB``  with the product of the Jacobian of the backward problem  at the
      point (``t``, ``y``, ``yB``) and the vector ``vB``.  Here, ``y`` is the
      solution of the original IVP at time ``t`` and  ``yB`` is the solution of
      the backward problem at the same time.  The rest of the arguments are
      equivalent to those passed to a function of type  ``IDALsJacTimesVecFn``
      (see :numref:`IDAS.Usage.SIM.user_supplied.jtimesFn`).

   .. warning::

      The previous function type ``IDASpilsJacTimesVecFnBS`` is identical to
      ``IDALsJacTimesVecFnBS``, and may still be used for
      backward-compatibility.  However, this will be deprecated in future
      releases, so we recommend that users transition to the new function type
      name soon.


.. _IDAS.Usage.ADJ.user_supplied.jactimesvecsetup_b:

Jacobian-vector product setup for the backward problem (matrix-free linear solvers)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If the user’s Jacobian-times-vector requires that any Jacobian-related data be
preprocessed or evaluated, then this needs to be done in a user-supplied
function of type :c:type:`IDALsJacTimesSetupFnB` or
:c:type:`IDALsJacTimesSetupFnBS`, defined as follows:

.. c:type:: int (*IDALsJacTimesSetupFnB)(realtype tt, N_Vector yy, N_Vector yp, N_Vector yB, N_Vector ypB, N_Vector resvalB, realtype cjB, void *user_dataB)

   This function preprocesses and/or evaluates Jacobian data needed  by the
   Jacobian-times-vector routine for the backward problem.

   **Arguments:**
     * ``tt`` -- is the current value of the independent variable.
     * ``yy`` -- is the current value of the dependent variable vector, :math:`y(t)` .
     * ``yp`` -- is the current value of :math:`\dot{y}(t)` .
     * ``yB`` -- is the current value of the backward dependent variable vector.
     * ``ypB`` -- is the current value of the backward dependent derivative vector.
     * ``resvalB`` -- is the current value of the residual for the backward problem.
     * ``cjB`` -- is the scalar in the system Jacobian, proportional to the inverse of the    step size ( :math:`\alpha` in :eq:`IDAS_DAE_Jacobian` ).
     * ``user_dataB`` -- is a pointer to user data — the same as the ``user_dataB`` parameter passed to :c:func:`IDASetUserDataB` .

   **Return value:**
      The value returned by the Jacobian-vector setup function
      should be  if successful, positive for a recoverable error (in
      which case the step will be retried), or negative for an
      unrecoverable error (in which case the integration is halted).

   **Notes:**
      Each call to the Jacobian-vector setup function is preceded by a call to
      the backward problem residual user function with the same  ``(t,y, yp, yB,
      ypB)`` arguments.  Thus, the setup function can use any auxiliary data
      that is computed  and saved during the evaluation of the DAE residual.  If
      the user's ``IDALsJacTimesVecFnB`` function uses difference quotient
      approximations, it may need to access quantities not in the call  list.
      These include the current stepsize, the error weights, etc.  To obtain
      these, the user will need to add a pointer to ``ida_mem``  to
      ``user_dataB`` and then use the ``IDAGet*`` functions described in
      :numref:`IDAS.Usage.SIM.user_callable.optional_output.main`.
      The unit roundoff can be accessed as  ``UNIT_ROUNDOFF`` defined in
      ``sundials_types.h``.

   .. warning::

      The previous function type ``IDASpilsJacTimesSetupFnB`` is identical to
      ``IDALsJacTimesSetupFnB``, and may still be used for
      backward-compatibility.  However, this will be deprecated in future
      releases, so we recommend that users transition to the new function type
      name soon.


.. c:type:: int (*IDALsJacTimesSetupFnBS)(realtype tt, N_Vector yy, N_Vector yp, N_Vector *yyS, N_Vector *ypS, N_Vector yB, N_Vector ypB, N_Vector resvalB, realtype cjB, void *user_dataB)

   This function preprocesses and/or evaluates Jacobian data needed  by the Jacobian-times-vector routine for the backward problem, in the case that  the backward problem depends on the forward sensitivities.

   **Arguments:**
     * ``tt`` -- is the current value of the independent variable.
     * ``yy`` -- is the current value of the dependent variable vector, :math:`y(t)` .
     * ``yp`` -- is the current value of :math:`\dot{y}(t)` .
     * ``yyS`` -- a pointer to an array of ``Ns`` vectors containing the sensitivities of    the forward solution.
     * ``ypS`` -- a pointer to an array of ``Ns`` vectors containing the derivatives of    the forward sensitivities.
     * ``yB`` -- is the current value of the backward dependent variable vector.
     * ``ypB`` -- is the current value of the backward dependent derivative vector.
     * ``resvalB`` -- is the current value of the residual for the backward problem.
     * ``cjB`` -- is the scalar in the system Jacobian, proportional to the inverse of the    step size ( :math:`\alpha` in :eq:`IDAS_DAE_Jacobian` ).
     * ``user_dataB`` -- is a pointer to user data — the same as the ``user_dataB`` parameter passed to :c:func:`IDASetUserDataB` .

   **Return value:**
      The value returned by the Jacobian-vector setup function should be  if
      successful, positive for a recoverable error (in which case the step will
      be retried), or negative for an unrecoverable error (in which case the
      integration is halted).

   **Notes:**
      Each call to the Jacobian-vector setup function is preceded by a call to
      the backward problem residual user function with the same  ``(t,y, yp,
      yyS, ypS, yB, ypB)`` arguments.  Thus, the setup function can use any
      auxiliary data that is computed  and saved during the evaluation of the
      DAE residual.  If the user's ``IDALsJacTimesVecFnB`` function uses
      difference quotient  approximations, it may need to access quantities not
      in the call  list. These include the current stepsize, the error weights,
      etc.  To obtain these, the user will need to add a pointer to ``ida_mem``
      to ``user_dataB`` and then use the ``IDAGet*`` functions described in
      :numref:`IDAS.Usage.ADJ.user_callable.optional_output_b`. The unit roundoff
      can be accessed as  ``UNIT_ROUNDOFF`` defined in ``sundials_types.h``.
      The previous function type ``IDASpilsJacTimesSetupFnBS`` is deprecated.

   .. warning::

      The previous function type ``IDASpilsJacTimesSetupFnBS`` is identical to
      :c:type:`IDALsJacTimesSetupFnBS`, and may still be used for
      backward-compatibility. However, this will be deprecated in future
      releases, so we recommend that users transition to the new function type
      name soon.


.. _IDAS.Usage.ADJ.user_supplied.psolve_b:

Preconditioner solve for the backward problem (iterative linear solvers)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If preconditioning is used during integration of the backward problem, then the
user must provide a function to solve the linear system :math:`Pz = r`, where
:math:`P` is a left preconditioner matrix. This function must have one of the
following two forms:

.. c:type:: int (*IDALsPrecSolveFnB)(realtype t, N_Vector yy, N_Vector yp, N_Vector yB, N_Vector ypB, N_Vector resvalB, N_Vector rvecB, N_Vector zvecB, realtype cjB, realtype deltaB, void *user_dataB)

   This function solves the preconditioning system :math:`Pz = r` for the backward
   problem.

   **Arguments:**
     * ``t`` -- is the current value of the independent variable.
     * ``yy`` -- is the current value of the forward solution vector.
     * ``yp`` -- is the current value of the forward solution derivative vector.
     * ``yB`` -- is the current value of the backward dependent variable vector.
     * ``ypB`` -- is the current value of the backward dependent derivative vector.
     * ``resvalB`` -- is the current value of the residual for the backward problem.
     * ``rvecB`` -- is the right-hand side vector :math:`r` of the linear system to be solved.
     * ``zvecB`` -- is the computed output vector.
     * ``cjB`` -- is the scalar in the system Jacobian, proportional to the inverse of the    step size ( :math:`\alpha` in :eq:`IDAS_DAE_Jacobian` ).
     * ``deltaB`` -- is an input tolerance to be used if an iterative method    is employed in the solution.
     * ``user_dataB`` -- is a pointer to user data — the same as the ``user_dataB`` parameter passed to the function :c:func:`IDASetUserDataB` .

   **Return value:**
      The return value of a preconditioner solve function for the backward
      problem should be  if successful, positive for a recoverable error (in
      which case the step will be retried), or negative for an unrecoverable
      error (in which case the integration is halted).

   .. warning::

      The previous function type ``IDASpilsPrecSolveFnB`` is identical to
      ``IDALsPrecSolveFnB``, and is deprecated.


.. c:type:: int (*IDALsPrecSolveFnBS)(realtype t, N_Vector yy, N_Vector yp, N_Vector *yyS, N_Vector *ypS, N_Vector yB, N_Vector ypB, N_Vector resvalB, N_Vector rvecB, N_Vector zvecB, realtype cjB, realtype deltaB, void *user_dataB)

   This function solves the preconditioning system :math:`Pz = r` for the
   backward problem,  for the case in which the backward problem depends on the
   forward sensitivities.

   **Arguments:**
     * ``t`` -- is the current value of the independent variable.
     * ``yy`` -- is the current value of the forward solution vector.
     * ``yp`` -- is the current value of the forward solution derivative vector.
     * ``yyS`` -- a pointer to an array of ``Ns`` vectors containing the sensitivities of    the forward solution.
     * ``ypS`` -- a pointer to an array of ``Ns`` vectors containing the derivatives of    the forward sensitivities.
     * ``yB`` -- is the current value of the backward dependent variable vector.
     * ``ypB`` -- is the current value of the backward dependent derivative vector.
     * ``resvalB`` -- is the current value of the residual for the backward problem.
     * ``rvecB`` -- is the right-hand side vector :math:`r` of the linear system to be solved.
     * ``zvecB`` -- is the computed output vector.
     * ``cjB`` -- is the scalar in the system Jacobian, proportional to the inverse of the    step size ( :math:`\alpha` in :eq:`IDAS_DAE_Jacobian` ).
     * ``deltaB`` -- is an input tolerance to be used if an iterative method    is employed in the solution.
     * ``user_dataB`` -- is a pointer to user data — the same as the ``user_dataB`` parameter passed to the function :c:func:`IDASetUserDataB` .

   **Return value:**
      The return value of a preconditioner solve function for the backward
      problem should be  if successful, positive for a recoverable error (in
      which case the step will be retried), or negative for an unrecoverable
      error (in which case the integration is halted).

   .. warning::

      The previous function type ``IDASpilsPrecSolveFnBS`` is identical to
      ``IDALsPrecSolveFnBS``, and is deprecated.


.. _IDAS.Usage.ADJ.user_supplied.psetup_b:

Preconditioner setup for the backward problem (iterative linear solvers)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If the user’s preconditioner requires that any Jacobian-related data be
preprocessed or evaluated, then this needs to be done in a user-supplied
function of one of the following two types:

.. c:type:: int (*IDALsPrecSetupFnB)(realtype t, N_Vector yy, N_Vector yp, N_Vector yB, N_Vector ypB, N_Vector resvalB, realtype cjB, void *user_dataB)

   This function preprocesses and/or evaluates Jacobian-related data needed  by
   the preconditioner for the backward problem.

   **Arguments:**
     * ``t`` -- is the current value of the independent variable.
     * ``yy`` -- is the current value of the forward solution vector.
     * ``yp`` -- is the current value of the forward solution vector.
     * ``yB`` -- is the current value of the backward dependent variable vector.
     * ``ypB`` -- is the current value of the backward dependent derivative vector.
     * ``resvalB`` -- is the current value of the residual for the backward problem.
     * ``cjB`` -- is the scalar in the system Jacobian, proportional to the inverse of the    step size ( :math:`\alpha` in :eq:`IDAS_DAE_Jacobian` ).
     * ``user_dataB`` -- is a pointer to user data — the same as the ``user_dataB`` parameter passed to the function :c:func:`IDASetUserDataB` .

   **Return value:**
      The return value of a preconditioner setup function for the backward
      problem should be  if successful, positive for a recoverable error (in
      which case the step will be retried), or negative for an unrecoverable
      error (in which case the integration is halted).

   .. warning::

      The previous function type ``IDASpilsPrecSetupFnB`` is identical to
      ``IDALsPrecSetupFnB``, and is deprecated.


.. c:type:: int (*IDALsPrecSetupFnBS)(realtype t, N_Vector yy, N_Vector yp, N_Vector *yyS, N_Vector *ypS, N_Vector yB, N_Vector ypB, N_Vector resvalB, realtype cjB, void *user_dataB)

   This function preprocesses and/or evaluates Jacobian-related data needed  by
   the preconditioner for the backward problem, in the case where the  backward
   problem depends on the forward sensitivities.

   **Arguments:**
     * ``t`` -- is the current value of the independent variable.
     * ``yy`` -- is the current value of the forward solution vector.
     * ``yp`` -- is the current value of the forward solution vector.
     * ``yyS`` -- a pointer to an array of ``Ns`` vectors containing the sensitivities of    the forward solution.
     * ``ypS`` -- a pointer to an array of ``Ns`` vectors containing the derivatives of    the forward sensitivities.
     * ``yB`` -- is the current value of the backward dependent variable vector.
     * ``ypB`` -- is the current value of the backward dependent derivative vector.
     * ``resvalB`` -- is the current value of the residual for the backward problem.
     * ``cjB`` -- is the scalar in the system Jacobian, proportional to the inverse of the    step size ( :math:`\alpha` in :eq:`IDAS_DAE_Jacobian` ).
     * ``user_dataB`` -- is a pointer to user data — the same as the ``user_dataB`` parameter passed to the function :c:func:`IDASetUserDataB` .

   **Return value:**
      The return value of a preconditioner setup function for the backward
      problem should be  if successful, positive for a recoverable error (in
      which case the step will be retried), or negative for an unrecoverable
      error (in which case the integration is halted).

   .. warning::

      The previous function type ``IDASpilsPrecSetupFnBS`` is identical to
      ``IDALsPrecSetupFnBS``, and is deprecated.


Using the band-block-diagonal preconditioner for backward problems
------------------------------------------------------------------

As on the forward integration phase, the efficiency of Krylov iterative methods
for the solution of linear systems can be greatly enhanced through
preconditioning. The band-block-diagonal preconditioner module IDABBDPRE,
provides interface functions through which it can be used on the backward
integration phase.

The adjoint module in IDAS offers an interface to the band-block-diagonal
preconditioner module IDABBDPRE described in section
:numref:`IDAS.Usage.precond.idabbdpre`. This generates a preconditioner that is
a block-diagonal matrix with each block being a band matrix and can be used with
one of the Krylov linear solvers and with the MPI-parallel vector module
``NVECTOR_PARALLEL``.

In order to use the IDABBDPRE module in the solution of the backward
problem, the user must define one or two additional functions, described at the
end of this section.

Usage of IDABBDPRE for the backward problem
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The IDABBDPRE module is initialized by calling the following function, *after*
an iterative linear solver for the backward problem has been attached to IDAS by
calling :c:func:`IDASetLinearSolverB` (see :numref:`IDAS.Usage.ADJ.user_callable.lin_solv_b`).

.. c:function:: int IDABBDPrecInitB(void * ida_mem, int which, sunindextype NlocalB, sunindextype mudqB, sunindextype mldqB, sunindextype mukeepB, sunindextype mlkeepB, realtype dqrelyB, IDABBDLocalFnB GresB, IDABBDCommFnB GcommB)

   The function :c:func:`IDABBDPrecInitB` initializes and allocates  memory for
   the IDABBDPRE preconditioner for the backward problem.

   **Arguments:**
     * ``ida_mem`` -- pointer to the IDAS memory block.
     * ``which`` -- the identifier of the backward problem.
     * ``NlocalB`` -- local vector dimension for the backward problem.
     * ``mudqB`` -- upper half-bandwidth to be used in the difference-quotient Jacobian approximation.
     * ``mldqB`` -- lower half-bandwidth to be used in the difference-quotient Jacobian approximation.
     * ``mukeepB`` -- upper half-bandwidth of the retained banded approximate Jacobian block.
     * ``mlkeepB`` -- lower half-bandwidth of the retained banded approximate Jacobian block.
     * ``dqrelyB`` -- the relative increment in components of ``yB`` used in the difference quotient
       approximations. The default is ``dqrelyB``  :math:`= \sqrt{\text{unit roundoff}}` , which can be
       specified by passing ``dqrely`` :math:`= 0.0`.
     * ``GresB`` -- the C function which computes :math:`G_B(t,y,\dot{y}, y_B, \dot{y}_B)`, the function approximating the residual of the backward problem.
     * ``GcommB`` -- the optional C function which performs all interprocess communication required for the computation of :math:`G_B`.

   **Return value:**
     * ``IDALS_SUCCESS`` -- The call to :c:func:`IDABBDPrecInitB` was successful.
     * ``IDALS_MEM_FAIL`` -- A memory allocation request has failed.
     * ``IDALS_MEM_NULL`` -- The ``ida_mem`` argument was ``NULL``.
     * ``IDALS_LMEM_NULL`` -- No linear solver has been attached.
     * ``IDALS_ILL_INPUT`` -- An invalid parameter has been passed.


To reinitialize the IDABBDPRE preconditioner module for the backward problem,
possibly with a change in ``mudqB``, ``mldqB``, or ``dqrelyB``, call the
following function:

.. c:function:: int IDABBDPrecReInitB(void * ida_mem, int which, sunindextype mudqB, sunindextype mldqB, realtype dqrelyB)

   The function :c:func:`IDABBDPrecReInitB` reinitializes the IDABBDPRE
   preconditioner  for the backward problem.

   **Arguments:**
     * ``ida_mem`` -- pointer to the IDAS memory block returned by :c:func:`IDACreate`.
     * ``which`` -- the identifier of the backward problem.
     * ``mudqB`` -- upper half-bandwidth to be used in the difference-quotient Jacobian approximation.
     * ``mldqB`` -- lower half-bandwidth to be used in the difference-quotient Jacobian approximation.
     * ``dqrelyB`` -- the relative increment in components of ``yB`` used in the difference quotient approximations.

   **Return value:**
     * ``IDALS_SUCCESS`` -- The call to :c:func:`IDABBDPrecReInitB` was successful.
     * ``IDALS_MEM_FAIL`` -- A memory allocation request has failed.
     * ``IDALS_MEM_NULL`` -- The ``ida_mem`` argument was ``NULL``.
     * ``IDALS_PMEM_NULL`` -- The :c:func:`IDABBDPrecInitB` has not been previously called.
     * ``IDALS_LMEM_NULL`` -- No linear solver has been attached.
     * ``IDALS_ILL_INPUT`` -- An invalid parameter has been passed.


User-supplied functions for IDABBDPRE
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

To use the IDABBDPRE module, the user must supply one or two functions which the
module calls to construct the preconditioner: a required function ``GresB`` (of
type :c:type:`IDABBDLocalFnB`) which approximates the residual of the backward
problem and which is computed locally, and an optional function ``GcommB`` (of
type :c:type:`IDABBDCommFnB`) which performs all interprocess communication
necessary to evaluate this approximate residual (see
:numref:`IDAS.Usage.precond.idabbdpre`).  The prototypes for these two functions
are described below.


.. c:type:: int (*IDABBDLocalFnB)(sunindextype NlocalB, realtype t, N_Vector y, N_Vector yp, N_Vector yB, N_Vector ypB, N_Vector gB, void *user_dataB)

   This ``GresB`` function loads the vector ``gB``, an approximation to the
   residual of the backward problem, as a function of ``t``, ``y``, ``yp``,  and
   ``yB`` and ``ypB``.

   **Arguments:**
     * ``NlocalB`` -- is the local vector length for the backward problem.
     * ``t`` -- is the value of the independent variable.
     * ``y`` -- is the current value of the forward solution vector.
     * ``yp`` -- is the current value of the forward solution derivative vector.
     * ``yB`` -- is the current value of the backward dependent variable vector.
     * ``ypB`` -- is the current value of the backward dependent derivative vector.
     * ``gB`` -- is the output vector, :math:`G_B(t,y,\dot y, y_B, \dot y_B)` .
     * ``user_dataB`` -- is a pointer to user data — the same as the ``user_dataB`` parameter passed to :c:func:`IDASetUserDataB` .

   **Return value:**
      An ``IDABBDLocalFnB`` should return 0 if successful, a positive value if a
      recoverable error occurred (in which case IDAS will attempt to correct),
      or a negative value if it failed unrecoverably (in which case the
      integration is halted and :c:func:`IDASolveB` returns
      ``IDA_LSETUP_FAIL``).

   **Notes:**
      This routine must assume that all interprocess communication of data
      needed to  calculate ``gB`` has already been done, and this data is
      accessible within  ``user_dataB``.

      .. warning::
        Before calling the user's ``IDABBDLocalFnB``, IDAS needs to evaluate
        (through interpolation) the values of the states from the forward
        integration.  If an error occurs in the interpolation, IDAS triggers an
        unrecoverable  failure in the preconditioner setup function which will
        halt the integration  (:c:func:`IDASolveB` returns ``IDA_LSETUP_FAIL``).


.. c:type:: int (*IDABBDCommFnB)(sunindextype NlocalB, realtype t, N_Vector y, N_Vector yp, N_Vector yB, N_Vector ypB, void *user_dataB)

   This ``GcommB`` function performs all interprocess communications necessary
   for the execution of the ``GresB`` function above, using the input  vectors
   ``y``, ``yp``, ``yB`` and ``ypB``.

   **Arguments:**
     * ``NlocalB`` -- is the local vector length.
     * ``t`` -- is the value of the independent variable.
     * ``y`` -- is the current value of the forward solution vector.
     * ``yp`` -- is the current value of the forward solution derivative vector.
     * ``yB`` -- is the current value of the backward dependent variable vector.
     * ``ypB`` -- is the current value of the backward dependent derivative vector.
     * ``user_dataB`` -- is a pointer to user data — the same as the ``user_dataB`` parameter passed to :c:func:`IDASetUserDataB` .

   **Return value:**
      An ``IDABBDCommFnB`` should return 0 if successful, a positive value if a
      recoverable error occurred (in which case IDAS will attempt to correct),
      or a negative value if it failed unrecoverably (in which case the
      integration is halted and :c:func:`IDASolveB` returns
      ``IDA_LSETUP_FAIL``).

   **Notes:**
      The ``GcommB`` function is expected to save communicated data in space
      defined within the  structure ``user_dataB``.

      Each call to the ``GcommB``
      function is preceded by a call to the function that  evaluates the
      residual of the backward problem with the same ``t``, ``y``, ``yp``,
      ``yB`` and ``ypB`` arguments. If there is no additional communication
      needed, then  pass ``GcommB`` :math:`=` ``NULL`` to
      :c:func:`IDABBDPrecInitB`.
