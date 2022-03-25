.. ----------------------------------------------------------------
   SUNDIALS Copyright Start
   Copyright (c) 2002-2022, Lawrence Livermore National Security
   and Southern Methodist University.
   All rights reserved.

   See the top-level LICENSE and NOTICE files for details.

   SPDX-License-Identifier: BSD-3-Clause
   SUNDIALS Copyright End
   ----------------------------------------------------------------

.. _CVODES.Usage.ADJ:

Using CVODES for Adjoint Sensitivity Analysis
=============================================

This chapter describes the use of CVODES to compute sensitivities of derived
functions using adjoint sensitivity analysis. As mentioned before, the adjoint
sensitivity module of CVODES provides the infrastructure for integrating
backward in time any system of ODEs that depends on the solution of the original
IVP, by providing various interfaces to the main CVODES integrator, as well as
several supporting user-callable functions. For this reason, in the following
sections we refer to the *backward problem* and not to the *adjoint problem*
when discussing details relevant to the ODEs that are integrated backward in
time. The backward problem can be the adjoint problem :eq:`CVODES_adj_eqns` or
:eq:`CVODES_adj_eqns`, and can be augmented with some quadrature differential
equations.

CVODES uses various constants for both input and output. These are defined as
needed in this chapter, but for convenience are also listed separately in
:numref:`CVODES.Constants`.

We begin with a brief overview, in the form of a skeleton user program.
Following that are detailed descriptions of the interface to the various
user-callable functions and of the user-supplied functions that were not already
described in :numref:`CVODES.Usage.SIM`.

.. _CVODES.Usage.ADJ.skeleton_sim:

A skeleton of the user’s main program
-------------------------------------

The following is a skeleton of the user’s main program as an application of
CVODES. The user program is to have these steps in the order indicated, unless
otherwise noted. For the sake of brevity, we defer many of the details to the
later sections. As in :numref:`CVODES.Usage.SIM.skeleton_sim`, most steps are
independent of the ``N_Vector``, ``SUNMatrix``, ``SUNLinearSolver``, and
``SUNNonlinearSolver`` implementations used. For the steps that are not, refer
to Chapters :numref:`NVectors`, :numref:`SUNMatrix`, :numref:`SUNLinSol`, and
:numref:`SUNNonlinSol` for the specific name of the function to be called or
macro to be referenced.

Steps that are unchanged from the skeleton programs presented in
:numref:`CVODES.Usage.SIM.skeleton_sim`, :numref:`CVODES.Usage.FSA.skeleton_sim`,
and :numref:`CVODES.Usage.purequad` are left unbolded.

#. Initialize parallel or multi-threaded environment, if appropriate

#. Create the SUNDIALS context object

#. Set problem dimensions etc. for the forward problem

#. Set initial conditions for the forward problem

#. Create CVODES object for the forward problem

#. Initialize CVODES for the forward problem

#. Specify integration tolerances for forward problem

#. Create matrix object for the forward problem

#. Create linear solver object for the forward problem

#. Set linear solver optional inputs for the forward problem

#. Attach linear solver module for the forward problem

#. Set optional inputs for the forward problem

#. Create nonlinear solver object for the forward problem

#. Attach nonlinear solver module for the forward problem

#. Set nonlinear solver optional inputs for the forward problem

#. Initialize quadrature problem or problems for forward problems, using :c:func:`CVodeQuadInit` and/or :c:func:`CVodeQuadSensInit`.

#. Initialize forward sensitivity problem

#. Specify rootfinding

#. **Allocate space for the adjoint computation**

   Call :c:func:`CVodeAdjInit` to allocate memory for the combined
   forward-backward problem. This call requires ``Nd``, the number of steps
   between two consecutive checkpoints. :c:func:`CVodeAdjInit` also specifies
   the type of interpolation used (see :numref:`CVODES.Mathematics.Checkpointing`).

#. **Integrate forward problem**

   Call :c:func:`CVodeF`, a wrapper for the CVODES main integration function
   :c:func:`CVode`, either in ``CV_NORMAL`` mode to the time ``tout`` or in
   ``CV_ONE_STEP`` mode inside a loop (if intermediate solutions of the forward
   problem are desired). The final value of ``tret`` is then the maximum
   allowable value for the endpoint :math:`T` of the backward problem.

#. **Set problem dimensions etc. for the backward problem**

   .. index:: back_start

   This generally includes the backward problem vector length ``NB``, and possibly the local vector length ``NBlocal``.

#. **Set initial values for the backward problem**

   Set the endpoint time ``tB0 = T``, and set the corresponding vector ``yB0``
   at which the backward problem starts.

#. **Create the backward problem**

   Call :c:func:`CVodeCreateB`, a wrapper for :c:func:`CVodeCreate`, to create
   the CVODES memory block for the new backward problem. Unlike
   :c:func:`CVodeCreate`, the function :c:func:`CVodeCreateB` does not return a
   pointer to the newly created memory block. Instead, this pointer is attached
   to the internal adjoint memory block (created by :c:func:`CVodeAdjInit`) and
   returns an identifier called ``which`` that the user must later specify in
   any actions on the newly created backward problem.

#. **Allocate memory for the backward problem**

   Call :c:func:`CVodeInitB` (or :c:func:`CVodeInitBS`, when the backward
   problem depends on the forward sensitivities). The two functions are actually
   wrappers for :c:func:`CVodeInit` and allocate internal memory, specify
   problem data, and initialize CVODES at ``tB0`` for the backward problem.

#. **Specify integration tolerances for backward problem**

   Call :c:func:`CVodeSStolerancesB` or :c:func:`CVodeSVtolerancesB` to specify
   a scalar relative tolerance and scalar absolute tolerance or scalar relative
   tolerance and a vector of absolute tolerances, respectively. The functions
   are wrappers for :c:func:`CVodeSStolerances` and :c:func:`CVodeSVtolerances`,
   but they require an extra argument ``which``, the identifier of the backward
   problem returned by :c:func:`CVodeCreateB`.

#. **Create matrix object for the backward problem**

   .. index:: matrixB

   If a nonlinear solver requiring a linear solve will be used (e.g., the the
   default Newton iteration) and the linear solver will be a direct linear
   solver, then a template Jacobian matrix must be created by calling the
   appropriate constructor function defined by the particular ``SUNMatrix``
   implementation.

   For the native SUNDIALS ``SUNMatrix`` implementations, the matrix object may
   be created using a call of the form ``SUN***Matrix(...)`` where ``***`` is
   the name of the matrix (see :numref:`SUNMatrix` for details).

#. **Create linear solver object for the backward problem**

   .. index:: lin_solverB

   If a nonlinear solver requiring a linear solver is chosen (e.g., the default
   Newton iteration), then the desired linear solver object for the backward
   problem must be created by calling the appropriate constructor function
   defined by the particular ``SUNLinearSolver`` implementation.

   For any of the SUNDIALS-supplied ``SUNLinearSolver`` implementations, the
   linear solver object may be created using a call of the form

   ``SUNLinearSolver LS = SUNLinSol_*(...);``

   where ``*`` can be replaced with “Dense”, “SPGMR”, or other options, as
   discussed in :numref:`CVODES.Usage.SIM.user_callable.lin_solv_init` and Chapter
   :numref:`SUNLinSol`.

   Note that it is not required to use the same linear solver module for both
   the forward and the backward problems; for example, the forward problem could
   be solved with the ``SUNLINSOL_BAND`` linear solver module and the backward
   problem with ``SUNLINSOL_SPGMR`` linear solver module.

#. **Set linear solver interface optional inputs for the backward problem**

   Call ``*Set*`` functions from the selected linear solver module to change
   optional inputs specific to that linear solver. See the documentation for
   each ``SUNLinearSolver`` module in Chapter :numref:`SUNLinSol`.

#. **Attach linear solver module for the backward problem**

   .. index:: lin_solver_interfaceB

   If a nonlinear solver requiring a linear solver is chosen for the backward
   problem (e.g., the default Newton iteration), then initialize the CVLS linear
   solver interface by attaching the linear solver object (and matrix object, if
   applicable) with the call to :c:func:`CVodeSetLinearSolverB`

   Alternately, if the CVODES-specific diagonal linear solver module, CVDIAG, is
   desired, initialize the linear solver module and attach it to CVODES with a
   call to :c:func:`CVDiagB`.

#. **Set optional inputs for the backward problem**

   Call ``CVodeSet*B`` functions to change from their default values any
   optional inputs that control the behavior of CVODES. Unlike their
   counterparts for the forward problem, these functions take an extra argument
   ``which``, the identifier of the backward problem returned by
   :c:func:`CVodeCreateB`.

#. **Create nonlinear solver object for the backward problem** (*optional*)

   If using a non-default nonlinear solver for the backward problem, then create
   the desired nonlinear solver object by calling the appropriate constructor
   function defined by the particular ``SUNNonlinearSolver`` implementation
   (e.g., ``NLSB = SUNNonlinSol_***(...);`` where ``***`` is the name of the
   nonlinear solver.

#. **Attach nonlinear solver module for the backward problem** (*optional*)

   If using a non-default nonlinear solver for the backward problem, then
   initialize the nonlinear solver interface by attaching the nonlinear
   solver object by calling :c:func:`CVodeSetNonlinearSolverB`.

#. **Initialize quadrature calculation**

   .. index:: quadB

   If additional quadrature equations must be evaluated, call
   :c:func:`CVodeQuadInitB` or :c:func:`CVodeQuadInitBS` (if quadrature depends
   also on the forward sensitivities). These functions are wrappers around
   :c:func:`CVodeQuadInit` and can be used to initialize and allocate memory for
   quadrature integration. Optionally, call ``CVodeSetQuad*B`` functions to
   change from their default values optional inputs that control the integration
   of quadratures during the backward phase.

#. **Integrate backward problem**

   Call :c:func:`CVodeB`, a second wrapper around the CVODES main integration
   function :c:func:`CVode`, to integrate the backward problem from ``tB0``.
   This function can be called either in ``CV_NORMAL`` or ``CV_ONE_STEP`` mode.
   Typically, :c:func:`CVodeB` will be called in ``CV_NORMAL`` mode with an end
   time equal to the initial time :math:`t_0` of the forward problem.

#. **Extract quadrature variables**

   .. index:: back_end

   If applicable, call :c:func:`CVodeGetQuadB`, a wrapper around
   :c:func:`CVodeGetQuad`, to extract the values of the quadrature variables at
   the time returned by the last call to :c:func:`CVodeB`.

#. **Deallocate memory**

   Upon completion of the backward integration, call all necessary deallocation
   functions. These include appropriate destructors for the vectors ``y`` and
   ``yB``, a call to :c:func:`CVodeFree` to free the CVODES memory block for the
   forward problem. If one or more additional Adjoint Sensitivity Analyses are
   to be done for this problem, a call to :c:func:`CVodeAdjFree` may be made to
   free and deallocate memory allocated for the backward problems, followed by a
   call to :c:func:`CVodeAdjInit`.

#. **Free the nonlinear solver memory for the forward and backward problems**

#. **Free linear solver and matrix memory for the forward and backward problems**

#. Free the SUNDIALS context with :c:func:`SUNContext_Free`

#. Finalize MPI, if used


The above user interface to the adjoint sensitivity module in CVODES was
motivated by the desire to keep it as close as possible in look and feel to the
one for ODE IVP integration. Note that if steps
:index:`back_start`-:index:`back_end` are not present, a program with the above
structure will have the same functionality as one described in
:numref:`CVODES.Usage.SIM.skeleton_sim` for integration of ODEs, albeit with some
overhead due to the checkpointing scheme.

If there are multiple backward problems associated with the same forward
problem, repeat steps :index:`back_start`-:index:`back_end` above for each
successive backward problem. In the process, each call to :c:func:`CVodeCreateB`
creates a new value of the identifier ``which``.


.. _CVODES.Usage.ADJ.user_callable:

User-callable functions for adjoint sensitivity analysis
--------------------------------------------------------


.. _CVODES.Usage.ADJ.user_callable.adjinit:

Adjoint sensitivity allocation and deallocation functions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

After the setup phase for the forward problem, but before the call to
:c:func:`CVodeF`, memory for the combined forward-backward problem must be
allocated by a call to the function :c:func:`CVodeAdjInit`. The form of the call
to this function is


.. c:function:: int CVodeAdjInit(void * cvode_mem, long int Nd, int interpType)

   The function :c:func:`CVodeAdjInit` updates CVODES memory block by allocating  the
   internal memory needed for backward integration.  Space is allocated for the
   ``Nd = N_d`` interpolation data points, and a linked  list of checkpoints
   is initialized.

   **Arguments:**
     * ``cvode_mem`` -- is the pointer to the CVODES memory block returned by a previous call to :c:func:`CVodeCreate`.
     * ``Nd`` -- is the number of integration steps between two consecutive checkpoints.
     * ``interpType`` -- specifies the type of interpolation used and can be ``CV_POLYNOMIAL`` or ``CV_HERMITE`` , indicating variable-degree polynomial and cubic Hermite interpolation, respectively see :numref:`CVODES.Mathematics.Checkpointing`.

   **Return value:**
     * ``CV_SUCCESS`` -- :c:func:`CVodeAdjInit` was successful.
     * ``CV_MEM_FAIL`` -- A memory allocation request has failed.
     * ``CV_MEM_NULL`` -- ``cvode_mem`` was NULL.
     * ``CV_ILL_INPUT`` -- One of the parameters was invalid: ``Nd`` was not positive or ``interpType`` is not one of the ``CV_POLYNOMIAL`` or ``CV_HERMITE``.

   **Notes:**
      The user must set ``Nd`` so that all data needed for interpolation of the
      forward problem solution between two checkpoints fits in memory.
      :c:func:`CVodeAdjInit`  attempts to allocate space for ``2*Nd+3`` variables
      of type ``N_Vector``.  If an error occurred, :c:func:`CVodeAdjInit` also sends a
      message to the  error handler function.


.. c:function:: int CVodeAdjReInit(void * cvode_mem)

   The function :c:func:`CVodeAdjReInit` reinitializes the CVODES memory  block for
   ASA, assuming that the number of steps between check  points and the type of
   interpolation remain unchanged.

   **Arguments:**
     * ``cvode_mem`` -- is the pointer to the CVODES memory block returned by a previous call to :c:func:`CVodeCreate`.

   **Return value:**
     * ``CV_SUCCESS`` -- :c:func:`CVodeAdjReInit` was successful.
     * ``CV_MEM_NULL`` -- ``cvode_mem`` was NULL.
     * ``CV_NO_ADJ`` -- The function :c:func:`CVodeAdjInit` was not previously called.

   **Notes:**
      The list of check points (and associated memory) is deleted.  The list of
      backward problems is kept. However, new backward problems can  be added to
      this list by calling :c:func:`CVodeCreateB`. If a new list of backward
      problems is also needed, then free the adjoint memory (by calling
      :c:func:`CVodeAdjFree`) and reinitialize ASA with :c:func:`CVodeAdjInit`.
      The CVODES memory for the forward and backward problems can be
      reinitialized  separately by calling :c:func:`CVodeReInit` and
      :c:func:`CVodeReInitB`, respectively.


.. c:function:: void CVodeAdjFree(void * cvode_mem)

   The function :c:func:`CVodeAdjFree` frees the memory related to backward
   integration allocated by a previous call to :c:func:`CVodeAdjInit`.

   **Argument:**
      * ``cvode_mem`` -- is the pointer to the CVODES memory block returned by a previous call to :c:func:`CVodeCreate`.

   **Return value:**
      The function has no return value.

   **Notes:**
      This function frees all memory allocated by :c:func:`CVodeAdjInit`. This
      includes workspace memory, the linked list of checkpoints, memory for the
      interpolation data, as well as the CVODES memory for the backward
      integration phase. Unless one or more further calls to
      :c:func:`CVodeAdjInit` are to be made, :c:func:`CVodeAdjFree` should not
      be called by the user, as it is invoked automatically by
      :c:func:`CVodeFree`.


.. _CVODES.Usage.ADJ.user_callable.cvodef:

Forward integration function
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The function :c:func:`CVodeF` is very similar to the CVODES function
:c:func:`CVode` in that it integrates the solution of the forward problem and
returns the solution in ``y``. At the same time, however, :c:func:`CVodeF`
stores checkpoint data every ``Nd`` integration steps. :c:func:`CVodeF` can be
called repeatedly by the user. Note that :c:func:`CVodeF` is used only for the
forward integration pass within an Adjoint Sensitivity Analysis. It is not for
use in Forward Sensitivity Analysis; for that, see :numref:`CVODES.Usage.FSA`.
The call to this function has the form


.. c:function:: int CVodeF(void * cvode_mem, realtype tout, N_Vector yret, realtype tret, int itask, int ncheck)

   The function :c:func:`CVodeF` integrates the forward problem over an interval
   in :math:`t`  and saves checkpointing data.

   **Arguments:**
     * ``cvode_mem`` -- pointer to the CVODES memory block.
     * ``tout`` -- the next time at which a computed solution is desired.
     * ``yret`` -- the computed solution vector :math:`y`.
     * ``tret`` -- the time reached by the solver output.
     * ``itask`` -- output mode a flag indicating the job of the solver for the next step. The ``CV_NORMAL`` task is to have the solver take internal steps until it has reached or just passed the user-specified ``tout`` parameter. The solver then interpolates in order to return an approximate value of :math:`y(\text{tout})`. The ``CV_ONE_STEP`` option tells the solver to just take one internal step and return the solution at the point reached by that step.
     * ``ncheck`` -- the number of internal checkpoints stored so far.

   **Return value:**
     * ``CV_SUCCESS`` -- :c:func:`CVodeF` succeeded.
     * ``CV_TSTOP_RETURN`` -- :c:func:`CVodeF` succeeded by reaching the optional stopping point.
     * ``CV_ROOT_RETURN`` -- :c:func:`CVodeF` succeeded and found one or more roots. In this case, ``tret`` is the location of the root. If ``nrtfn > 1`` , call :c:func:`CVodeGetRootInfo` to see which :math:`g_i` were found to have a root.
     * ``CV_NO_MALLOC`` -- The function :c:func:`CVodeInit` has not been previously called.
     * ``CV_ILL_INPUT`` -- One of the inputs to :c:func:`CVodeF` is illegal.
     * ``CV_TOO_MUCH_WORK`` -- The solver took ``mxstep`` internal steps but could not reach ``tout``.
     * ``CV_TOO_MUCH_ACC`` -- The solver could not satisfy the accuracy demanded by the user for some internal step.
     * ``CV_ERR_FAILURE`` -- Error test failures occurred too many times during one internal time step or occurred with :math:`|h| = h_{min}`.
     * ``CV_CONV_FAILURE`` -- Convergence test failures occurred too many times during one internal time step or occurred with :math:`|h| = h_{min}`.
     * ``CV_LSETUP_FAIL`` -- The linear solver's setup function failed in an unrecoverable manner.
     * ``CV_LSOLVE_FAIL`` -- The linear solver's solve function failed in an unrecoverable manner.
     * ``CV_NO_ADJ`` -- The function :c:func:`CVodeAdjInit` has not been previously called.
     * ``CV_MEM_FAIL`` -- A memory allocation request has failed in an attempt to allocate space for a new checkpoint.

   **Notes:**
      All failure return values are negative and therefore a test
      ``flag``:math:`< 0`  will trap all :c:func:`CVodeF` failures.  At this
      time, :c:func:`CVodeF` stores checkpoint information in memory only.
      Future versions will provide for a safeguard option of dumping checkpoint
      data into a temporary file as needed. The data stored at each checkpoint
      is basically  a snapshot of the CVODES internal memory block and contains
      enough information  to restart the integration from that time and to
      proceed with the same step size and  method order sequence as during the
      forward integration.  In addition, :c:func:`CVodeF` also stores
      interpolation data between consecutive checkpoints  so that, at the end of
      this first forward integration phase, interpolation information  is
      already available from the last checkpoint forward. In particular,  if no
      checkpoints were necessary, there is no need for the second forward
      integration phase.

   .. warning::
      It is illegal to change the integration tolerances between consecutive
      calls  to :c:func:`CVodeF`, as this information is not captured in the
      checkpoint data.


.. _CVODES.Usage.ADJ.user_callable.cvinitb:

Backward problem initialization functions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The functions :c:func:`CVodeCreateB` and :c:func:`CVodeInitB` (or
:c:func:`CVodeInitBS`) must be called in the order listed. They instantiate a
CVODES solver object, provide problem and solution specifications, and allocate
internal memory for the backward problem.


.. c:function:: int CVodeCreateB(void * cvode_mem, int lmmB, int which)

   The function :c:func:`CVodeCreateB` instantiates a CVODES solver object and
   specifies  the solution method for the backward problem.

   **Arguments:**
     * ``cvode_mem`` -- pointer to the CVODES memory block returned by :c:func:`CVodeCreate`.
     * ``lmmB`` -- specifies the linear multistep method and may be one of two possible values: ``CV_ADAMS`` or ``CV_BDF``.
     * ``which`` -- contains the identifier assigned by CVODES for the newly created backward problem. Any call to ``CVode*B`` functions requires such an identifier.

   **Return value:**
     * ``CV_SUCCESS`` -- The call to :c:func:`CVodeCreateB` was successful.
     * ``CV_MEM_NULL`` -- ``cvode_mem`` was ``NULL``.
     * ``CV_NO_ADJ`` -- The function :c:func:`CVodeAdjInit` has not been previously called.
     * ``CV_MEM_FAIL`` -- A memory allocation request has failed.


There are two initialization functions for the backward problem – one for the
case when the backward problem does not depend on the forward sensitivities, and
one for the case when it does. These two functions are described next.


.. c:function:: int CVodeInitB(void * cvode_mem, int which, CVRhsFnB rhsB, realtype tB0, N_Vector yB0)

   The function :c:func:`CVodeInitB` provides problem specification, allocates
   internal memory,  and initializes the backward problem.

   **Arguments:**
     * ``cvode_mem`` -- pointer to the CVODES memory block returned by :c:func:`CVodeCreate`.
     * ``which`` -- represents the identifier of the backward problem.
     * ``rhsB`` -- is the :c:type:`CVRhsFnB` function which computes :math:`f_B` , the right-hand side of the backward ODE problem.
     * ``tB0`` -- specifies the endpoint :math:`T` where final conditions are provided for the backward problem, normally equal to the endpoint of the forward integration.
     * ``yB0`` -- is the initial value at :math:`t =` ``tB0`` of the backward solution.

   **Return value:**
     * ``CV_SUCCESS`` -- The call to :c:func:`CVodeInitB` was successful.
     * ``CV_NO_MALLOC`` -- The function :c:func:`CVodeInit` has not been previously called.
     * ``CV_MEM_NULL`` -- ``cvode_mem`` was ``NULL``.
     * ``CV_NO_ADJ`` -- The function :c:func:`CVodeAdjInit` has not been previously called.
     * ``CV_BAD_TB0`` -- The final time ``tB0`` was outside the interval over which the forward problem was solved.
     * ``CV_ILL_INPUT`` -- The parameter ``which`` represented an invalid identifier, or either ``yB0`` or ``rhsB`` was ``NULL``.

   **Notes:**
      The memory allocated by :c:func:`CVodeInitB` is deallocated by the
      function  :c:func:`CVodeAdjFree`.


The function :c:func:`CVodeInitB` initializes the backward problem when it does
not depend on the forward sensitivities. It is essentially a wrapper for
:c:func:`CVodeInit` with some particularization for backward integration, as
described below.

For the case when backward problem also depends on the forward sensitivities,
user must call :c:func:`CVodeInitBS` instead of :c:func:`CVodeInitB`. Only the
third argument of each function differs between these two functions.


.. c:function:: int CVodeInitBS(void * cvode_mem, int which, CVRhsFnBS rhsBS, realtype tB0, N_Vector yB0)

   The function :c:func:`CVodeInitBS` provides problem specification, allocates
   internal memory,  and initializes the backward problem.

   **Arguments:**
     * ``cvode_mem`` -- pointer to the CVODES memory block returned by :c:func:`CVodeCreate`.
     * ``which`` -- represents the identifier of the backward problem.
     * ``rhsBS`` -- is the :c:type:`CVRhsFnBS` function which computes :math:`f_B` , the right-hand side of the backward ODE problem.
     * ``tB0`` -- specifies the endpoint :math:`T` where final conditions are provided for the backward problem.
     * ``yB0`` -- is the initial value at :math:`t =` ``tB0`` of the backward solution.

   **Return value:**
     * ``CV_SUCCESS`` -- The call to :c:func:`CVodeInitB` was successful.
     * ``CV_NO_MALLOC`` -- The function :c:func:`CVodeInit` has not been previously called.
     * ``CV_MEM_NULL`` -- ``cvode_mem`` was ``NULL``.
     * ``CV_NO_ADJ`` -- The function :c:func:`CVodeAdjInit` has not been previously called.
     * ``CV_BAD_TB0`` -- The final time ``tB0`` was outside the interval over which the forward problem was solved.
     * ``CV_ILL_INPUT`` -- The parameter ``which`` represented an invalid identifier, either ``yB0`` or ``rhsBS`` was ``NULL`` , or sensitivities were not active during the forward integration.

   **Notes:**
      The memory allocated by :c:func:`CVodeInitBS` is deallocated by the
      function  :c:func:`CVodeAdjFree`.


The function :c:func:`CVodeReInitB` reinitializes CVODES for the solution of a
series of backward problems, each identified by a value of the parameter
``which``. :c:func:`CVodeReInitB` is essentially a wrapper for
:c:func:`CVodeReInit`, and so all details given for :c:func:`CVodeReInit` apply
here. Also note that :c:func:`CVodeReInitB` can be called to reinitialize the
backward problem even it has been initialized with the sensitivity-dependent
version :c:func:`CVodeInitBS`. Before calling :c:func:`CVodeReInitB` for a new
backward problem, call any desired solution extraction functions ``CVodeGet**``
associated with the previous backward problem. The call to the
:c:func:`CVodeReInitB` function has the form


.. c:function:: int CVodeReInitB(void * cvode_mem, int which, realtype tB0, N_Vector yB0)

   The function :c:func:`CVodeReInitB` reinitializes a CVODES backward problem.

   **Arguments:**
     * ``cvode_mem`` -- pointer to CVODES memory block returned by :c:func:`CVodeCreate`.
     * ``which`` -- represents the identifier of the backward problem.
     * ``tB0`` -- specifies the endpoint :math:`T` where final conditions are provided for the backward problem.
     * ``yB0`` -- is the initial value at :math:`t =` ``tB0`` of the backward solution.

   **Return value:**
     * ``CV_SUCCESS`` -- The call to :c:func:`CVodeReInitB` was successful.
     * ``CV_NO_MALLOC`` -- The function :c:func:`CVodeInit` has not been previously called.
     * ``CV_MEM_NULL`` -- The ``cvode_mem`` memory block pointer was ``NULL``.
     * ``CV_NO_ADJ`` -- The function :c:func:`CVodeAdjInit` has not been previously called.
     * ``CV_BAD_TB0`` -- The final time ``tB0`` is outside the interval over which the forward problem was solved.
     * ``CV_ILL_INPUT`` -- The parameter ``which`` represented an invalid identifier, or ``yB0`` was ``NULL``.


.. _CVODES.Usage.ADJ.user_callable.cvtolerances_b:

Tolerance specification functions for backward problem
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

One of the following two functions must be called to specify the integration
tolerances for the backward problem. Note that this call must be made after the
call to :c:func:`CVodeInitB` or :c:func:`CVodeInitBS`.


.. c:function:: int CVodeSStolerancesB(void * cvode_mem, int which, realtype reltolB, realtype abstolB)

   The function :c:func:`CVodeSStolerancesB` specifies scalar relative and absolute  tolerances.

   **Arguments:**
     * ``cvode_mem`` -- pointer to the CVODES memory block returned by :c:func:`CVodeCreate`.
     * ``which`` -- represents the identifier of the backward problem.
     * ``reltolB`` -- is the scalar relative error tolerance.
     * ``abstolB`` -- is the scalar absolute error tolerance.

   **Return value:**
     * ``CV_SUCCESS`` -- The call to :c:func:`CVodeSStolerancesB` was successful.
     * ``CV_MEM_NULL`` -- The CVODES memory block was not initialized through a previous call to :c:func:`CVodeCreate`.
     * ``CV_NO_MALLOC`` -- The allocation function :c:func:`CVodeInit` has not been called.
     * ``CV_NO_ADJ`` -- The function :c:func:`CVodeAdjInit` has not been previously called.
     * ``CV_ILL_INPUT`` -- One of the input tolerances was negative.


.. c:function:: int CVodeSVtolerancesB(void * cvode_mem, int which, reltolBabstolB)

   The function :c:func:`CVodeSVtolerancesB` specifies scalar relative tolerance and  vector absolute tolerances.

   **Arguments:**
     * ``cvode_mem`` -- pointer to the CVODES memory block returned by :c:func:`CVodeCreate`.
     * ``which`` -- represents the identifier of the backward problem.
     * ``reltol`` -- is the scalar relative error tolerance.
     * ``abstol`` -- is the vector of absolute error tolerances.

   **Return value:**
     * ``CV_SUCCESS`` -- The call to :c:func:`CVodeSVtolerancesB` was successful.
     * ``CV_MEM_NULL`` -- The CVODES memory block was not initialized through a previous call to :c:func:`CVodeCreate`.
     * ``CV_NO_MALLOC`` -- The allocation function :c:func:`CVodeInit` has not been called.
     * ``CV_NO_ADJ`` -- The function :c:func:`CVodeAdjInit` has not been previously called.
     * ``CV_ILL_INPUT`` -- The relative error tolerance was negative or the absolute tolerance had a negative component.

   **Notes:**
      This choice of tolerances is important when the absolute error tolerance
      needs to  be different for each component of the state vector :math:`y`.


.. _CVODES.Usage.ADJ.user_callable.lin_solv_b:

Linear solver initialization functions for backward problem
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

All CVODES linear solver modules available for forward problems are available
for the backward problem. They should be created as for the forward problem and
then attached to the memory structure for the backward problem using the
following functions.


.. c:function:: int CVodeSetLinearSolverB(void * cvode_mem, int which, SUNLinearSolver LS, SUNMatrix A)

   The function :c:func:`CVodeSetLinearSolverB` attaches a generic
   ``SUNLinearSolver`` object ``LS`` and corresponding template Jacobian
   ``SUNMatrix`` object ``A`` to CVODES, initializing the  CVLS linear solver
   interface for solution of the backward  problem.

   **Arguments:**
     * ``cvode_mem`` -- pointer to the CVODES memory block.
     * ``which`` -- represents the identifier of the backward problem returned by :c:func:`CVodeCreateB`.
     * ``LS`` -- SUNLINSOL object to use for solving linear systems for the backward problem.
     * ``A`` -- SUNMATRIX object for used as a template for the Jacobian for the backward problem or ``NULL`` if not applicable.

   **Return value:**
     * ``CVLS_SUCCESS`` -- The CVLS initialization was successful.
     * ``CVLS_MEM_NULL`` -- The ``cvode_mem`` pointer is ``NULL``.
     * ``CVLS_ILL_INPUT`` -- The parameter ``which`` represented an invalid identifier.
     * ``CVLS_MEM_FAIL`` -- A memory allocation request failed.
     * ``CVLS_NO_ADJ`` -- The function ``CVAdjInit`` has not been previously called.

   **Notes:**
      If ``LS`` is a matrix-based linear solver, then the template  Jacobian
      matrix ``J`` will be used in the solve process, so if  additional storage
      is required within the ``SUNMatrix`` object  (e.g., for factorization of a
      banded matrix), ensure that the input  object is allocated with sufficient
      size (see the documentation of  the particular ``SUNMatrix`` type in
      :numref:`SUNMatrix`).  The previous routines ``CVDlsSetLinearSolverB`` and
      ``CVSpilsSetLinearSolverB`` are now wrappers for this routine, and may
      still be used for backward-compatibility.  However, these will be
      deprecated in future releases, so we recommend that users transition  to
      the new routine name soon.


.. c:function:: int CVDiagB(void * cvode_mem, int which)

   The function ``CVDiagB`` selects the CVDIAG linear solver for the solution
   of the backward problem.  The user's main program must include the
   ``cvodes_diag.h`` header file.

   **Arguments:**
     * ``cvode_mem`` -- pointer to the CVODES memory block.
     * ``which`` -- represents the identifier of the backward problem returned by :c:func:`CVodeCreateB`.

   **Return value:**
     * ``CVDIAG_SUCCESS`` -- The CVDIAG initialization was successful.
     * ``CVDIAG_MEM_NULL`` -- The ``cvode_mem`` pointer is ``NULL``.
     * ``CVDIAG_ILL_INPUT`` -- The CVDIAG solver is not compatible with the current NVECTOR module.
     * ``CVDIAG_MEM_FAIL`` -- A memory allocation request failed.

   **Notes:**
      The CVDIAG solver is the simplest of all of the available CVODES  linear
      solver interfaces.  The CVDIAG solver uses an approximate  diagonal
      Jacobian formed by way of a difference quotient. The user  does not have
      the option of supplying a function to compute an  approximate diagonal
      Jacobian.


.. _CVODES.Usage.ADJ.user_callable.nonlin_solv_init_b:

Nonlinear solver initialization function for backward problem
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

All CVODES nonlinear solver modules available for forward problems are available
for the backward problem. As with the forward problem CVODES uses the
``SUNNonlinearSolver`` implementation of Newton’s method defined by the
:ref:`SUNNONLINSOL_NEWTON <SUNNonlinSol.Newton>` module by default.

To specify a different nonlinear solver for the backward problem, the user’s
program must create a ``SUNNonlinearSolver`` object by calling the appropriate
constructor routine. The user must then attach the ``SUNNonlinearSolver`` object
by calling :c:func:`CVodeSetNonlinearSolverB`, as documented below.

When changing the nonlinear solver in CVODES, :c:func:`CVodeSetNonlinearSolverB`
must be called after :c:func:`CVodeInitB`. If any calls to :c:func:`CVodeB` have
been made, then CVODES will need to be reinitialized by calling
:c:func:`CVodeReInitB` to ensure that the nonlinear solver is initialized
correctly before any subsequent calls to :c:func:`CVodeB`.

.. c:function:: int CVodeSetNonlinearSolverB(void * cvode_mem, int which, SUNNonlinearSolver NLS)

   The function :c:func:`CVodeSetNonLinearSolverB` attaches a
   ``SUNNONLINEARSOLVER``  object (``NLS``) to CVODES for the solution of the
   backward problem.

   **Arguments:**
     * ``cvode_mem`` -- pointer to the CVODES memory block.
     * ``which`` -- represents the identifier of the backward problem returned by :c:func:`CVodeCreateB`.
     * ``NLS`` -- SUNNONLINSOL object to use for solving nonlinear systems for the backward problem.

   **Return value:**
     * ``CV_SUCCESS`` -- The nonlinear solver was successfully attached.
     * ``CV_MEM_NULL`` -- The ``cvode_mem`` pointer is ``NULL``.
     * ``CVLS_NO_ADJ`` -- The function ``CVAdjInit`` has not been previously called.
     * ``CV_ILL_INPUT`` -- The parameter ``which`` represented an invalid identifier or the SUNNONLINSOL object is ``NULL`` , does not implement the required nonlinear solver operations, is not of the correct type, or the residual function, convergence test function, or maximum number of nonlinear iterations could not be set.


.. _CVODES.Usage.ADJ.user_callable.cvsolveb:

Backward integration function
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The function :c:func:`CVodeB` performs the integration of the backward problem.
It is essentially a wrapper for the CVODES main integration function
:c:func:`CVode` and, in the case in which checkpoints were needed, it evolves
the solution of the backward problem through a sequence of forward-backward
integration pairs between consecutive checkpoints. The first run of each pair
integrates the original IVP forward in time and stores interpolation data; the
second run integrates the backward problem backward in time and performs the
required interpolation to provide the solution of the IVP to the backward
problem.

The function :c:func:`CVodeB` does not return the solution ``yB`` itself. To
obtain that, call the function :c:func:`CVodeGetB`, which is also described
below.

The :c:func:`CVodeB` function does not support rootfinding, unlike
:c:func:`CVodeF`, which supports the finding of roots of functions of
:math:`(t,y)`. If rootfinding was performed by :c:func:`CVodeF`, then for the
sake of efficiency, it should be disabled for :c:func:`CVodeB` by first calling
:c:func:`CVodeRootInit` with ``nrtfn`` = 0.

The call to :c:func:`CVodeB` has the form

.. c:function:: int CVodeB(void * cvode_mem, realtype tBout, int itaskB)

   The function :c:func:`CVodeB` integrates the backward ODE problem.

   **Arguments:**
     * ``cvode_mem`` -- pointer to the CVODES memory returned by :c:func:`CVodeCreate`.
     * ``tBout`` -- the next time at which a computed solution is desired.
     * ``itaskB`` --  output mode a flag indicating the job of the solver for the next step. The ``CV_NORMAL`` task is to have the solver take internal steps until it has reached or just passed the user-specified value ``tBout``. The solver then interpolates in order to return an approximate value of :math:`yB(\texttt{tBout})`. The ``CV_ONE_STEP`` option tells the solver to take just one internal step in the direction of ``tBout`` and return.

   **Return value:**
     * ``CV_SUCCESS`` -- :c:func:`CVodeB` succeeded.
     * ``CV_MEM_NULL`` -- ``cvode_mem`` was ``NULL``.
     * ``CV_NO_ADJ`` -- The function :c:func:`CVodeAdjInit` has not been previously called.
     * ``CV_NO_BCK`` -- No backward problem has been added to the list of backward problems by a call to :c:func:`CVodeCreateB`.
     * ``CV_NO_FWD`` -- The function :c:func:`CVodeF` has not been previously called.
     * ``CV_ILL_INPUT`` -- One of the inputs to :c:func:`CVodeB` is illegal.
     * ``CV_BAD_ITASK`` -- The ``itaskB`` argument has an illegal value.
     * ``CV_TOO_MUCH_WORK`` -- The solver took ``mxstep`` internal steps but could not reach ``tBout``.
     * ``CV_TOO_MUCH_ACC`` -- The solver could not satisfy the accuracy demanded by the user for some internal step.
     * ``CV_ERR_FAILURE`` -- Error test failures occurred too many times during one internal time step.
     * ``CV_CONV_FAILURE`` -- Convergence test failures occurred too many times during one internal time step.
     * ``CV_LSETUP_FAIL`` -- The linear solver's setup function failed in an unrecoverable manner.
     * ``CV_SOLVE_FAIL`` -- The linear solver's solve function failed in an unrecoverable manner.
     * ``CV_BCKMEM_NULL`` -- The solver memory for the backward problem was not created with a call to :c:func:`CVodeCreateB`.
     * ``CV_BAD_TBOUT`` -- The desired output time ``tBout`` is outside the interval over which the forward problem was solved.
     * ``CV_REIFWD_FAIL`` -- Reinitialization of the forward problem failed at the first checkpoint corresponding to the initial time of the forward problem.
     * ``CV_FWD_FAIL`` -- An error occurred during the integration of the forward problem.

   **Notes:**
      All failure return values are negative and therefore a test
      ``flag < 0``  will trap all :c:func:`CVodeB` failures.  In the case
      of multiple checkpoints and multiple backward problems, a given  call to
      :c:func:`CVodeB` in ``CV_ONE_STEP`` mode may not advance every problem
      one step, depending on the relative locations of the current times
      reached.  But repeated calls will eventually advance all problems to
      ``tBout``.


In the case of multiple checkpoints and multiple backward problems, a given call
to :c:func:`CVodeB` in ``CV_ONE_STEP`` mode may not advance every problem one
step, depending on the relative locations of the current times reached. But
repeated calls will eventually advance all problems to ``tBout``.

To obtain the solution ``yB`` to the backward problem, call the function :c:func:`CVodeGetB` as follows:

.. c:function:: int CVodeGetB(void * cvode_mem, int which, realtype tret, N_Vector yB)

   The function :c:func:`CVodeGetB` provides the solution ``yB`` of the backward
   ODE  problem.

   **Arguments:**
     * ``cvode_mem`` -- pointer to the CVODES memory returned by :c:func:`CVodeCreate`.
     * ``which`` -- the identifier of the backward problem.
     * ``tret`` -- the time reached by the solver output.
     * ``yB`` -- the backward solution at time ``tret``.

   **Return value:**
     * ``CV_SUCCESS`` -- :c:func:`CVodeGetB` was successful.
     * ``CV_MEM_NULL`` -- ``cvode_mem`` is ``NULL``.
     * ``CV_NO_ADJ`` -- The function :c:func:`CVodeAdjInit` has not been previously called.
     * ``CV_ILL_INPUT`` -- The parameter ``which`` is an invalid identifier.

   .. warning::
      The user must allocate space for ``yB``.  To obtain the solution
      associated with a given backward problem at some  other time within the
      last integration step, first obtain a pointer to the  proper CVODES
      memory structure by calling :c:func:`CVodeGetAdjCVodeBmem`  and then
      use it to call :c:func:`CVodeGetDky`.


Adjoint sensitivity optional input
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

At any time during the integration of the forward problem, the user can disable
the checkpointing of the forward sensitivities by calling the following
function:

.. c:function:: int CVodeAdjSetNoSensi(void * cvode_mem)

   The function :c:func:`CVodeAdjSetNoSensi` instructs :c:func:`CVodeF` not  to
   save checkpointing data for forward sensitivities anymore.

   **Arguments:**
     * ``cvode_mem`` -- pointer to the CVODES memory block.

   **Return value:**
     * ``CV_SUCCESS`` -- The call to :c:func:`CVodeCreateB` was successful.
     * ``CV_MEM_NULL`` -- ``cvode_mem`` was ``NULL``.
     * ``CV_NO_ADJ`` -- The function :c:func:`CVodeAdjInit` has not been previously called.


.. _CVODES.Usage.ADJ.user_callable.optional_input_b:

Optional input functions for the backward problem
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

As for the forward problem there are numerous optional input parameters that
control the behavior of the CVODES solver for the backward problem. CVODES
provides functions that can be used to change these optional input parameters
from their default values which are then described in detail in the remainder of
this section, beginning with those for the main CVODES solver and continuing
with those for the linear solver interfaces. Note that the diagonal linear
solver module has no optional inputs. For the most casual use of CVODES, the
reader can skip to :numref:`CVODES.Usage.ADJ.user_supplied`.

We note that, on an error return, all of the optional input functions send an
error message to the error handler function. All error return values are
negative, so the test ``flag < 0`` will catch all errors. Finally, a call to a
``CVodeSet***B`` function can be made from the user’s calling program at any
time and, if successful, takes effect immediately.

Main solver optional input functions
""""""""""""""""""""""""""""""""""""

The adjoint module in CVODES provides wrappers for most of the optional input
functions defined in :numref:`CVODES.Usage.SIM.optional_input.optin_main`. The
only difference is that the user must specify the identifier ``which`` of the
backward problem within the list managed by CVODES.

The optional input functions defined for the backward problem are:

.. code-block:: c

     flag = CVodeSetUserDataB(cvode_mem, which, user_dataB);
     flag = CVodeSetMaxOrdB(cvode_mem, which, maxordB);
     flag = CVodeSetMaxNumStepsB(cvode_mem, which, mxstepsB);
     flag = CVodeSetInitStepB(cvode_mem, which, hinB)
     flag = CVodeSetMinStepB(cvode_mem, which, hminB);
     flag = CVodeSetMaxStepB(cvode_mem, which, hmaxB);
     flag = CVodeSetStabLimDetB(cvode_mem, which, stldetB);
     flag = CVodeSetConstraintsB(cvode_mem, which, constraintsB);

Their return value ``flag`` (of type ``int``) can have any of the return values
of their counterparts, but it can also be ``CV_NO_ADJ`` if
:c:func:`CVodeAdjInit` has not been called, or ``CV_ILL_INPUT`` if ``which`` was
an invalid identifier.

Linear solver interface optional input functions
""""""""""""""""""""""""""""""""""""""""""""""""

When using matrix-based linear solver modules, the CVLS solver interface needs a
function to compute an approximation to the Jacobian matrix or the linear system
for the backward problem. The function to evaluate the Jacobian can be attached
through a call to either :c:func:`CVodeSetJacFnB` or :c:func:`CVodeSetJacFnBS`,
with the second used when the backward problem depends on the forward
sensitivities.

.. c:function:: int CVodeSetJacFnB(void * cvode_mem, int which, CVLsJacFnB jacB)

   The function :c:func:`CVodeSetJacFnB` specifies the Jacobian  approximation
   function to be used for the backward problem.

   **Arguments:**
     * ``cvode_mem`` -- pointer to the CVODES memory returned by :c:func:`CVodeCreate`.
     * ``which`` -- represents the identifier of the backward problem.
     * ``jacB`` -- user-defined Jacobian approximation function.

   **Return value:**
     * ``CVLS_SUCCESS`` -- :c:func:`CVodeSetJacFnB` succeeded.
     * ``CVLS_MEM_NULL`` -- ``cvode_mem`` was ``NULL``.
     * ``CVLS_NO_ADJ`` -- The function :c:func:`CVodeAdjInit` has not been previously called.
     * ``CVLS_LMEM_NULL`` -- The linear solver has not been initialized with a call to :c:func:`CVodeSetLinearSolverB`.
     * ``CVLS_ILL_INPUT`` -- The parameter ``which`` represented an invalid identifier.

   **Notes:**
      The previous routine :c:type:`CVDlsSetJacFnB` is now deprecated.


.. c:function:: int CVodeSetJacFnBS(void * cvode_mem, int which, CVLsJacFnBS jacBS)

   The function :c:func:`CVodeSetJacFnBS` specifies the Jacobian  approximation
   function to be used for the backward problem, in the  case where the backward
   problem depends on the forward sensitivities.

   **Arguments:**
     * ``cvode_mem`` -- pointer to the CVODES memory returned by :c:func:`CVodeCreate`.
     * ``which`` -- represents the identifier of the backward problem.
     * ``jacBS`` -- user-defined Jacobian approximation function.

   **Return value:**
     * ``CVLS_SUCCESS`` -- :c:func:`CVodeSetJacFnBS` succeeded.
     * ``CVLS_MEM_NULL`` -- ``cvode_mem`` was ``NULL``.
     * ``CVLS_NO_ADJ`` -- The function :c:func:`CVodeAdjInit` has not been previously called.
     * ``CVLS_LMEM_NULL`` -- The linear solver has not been initialized with a call to :c:func:`CVodeSetLinearSolverB`.
     * ``CVLS_ILL_INPUT`` -- The parameter ``which`` represented an invalid identifier.

   **Notes:**
      The previous routine :c:type:`CVDlsSetJacFnBS` is now deprecated.


.. c:function:: int CVodeSetLinSysFnB(void * cvode_mem, int which, CVLsLinSysFnB linsysB)

   The function :c:func:`CVodeSetLinSysFnB` specifies the linear system
   approximation function to be used for the backward problem.

   **Arguments:**
     * ``cvode_mem`` -- pointer to the CVODES memory returned by :c:func:`CVodeCreate`.
     * ``which`` -- represents the identifier of the backward problem.
     * ``linsysB`` -- user-defined linear system approximation function.

   **Return value:**
     * ``CVLS_SUCCESS`` -- :c:func:`CVodeSetLinSysFnB` succeeded.
     * ``CVLS_MEM_NULL`` -- ``cvode_mem`` was ``NULL``.
     * ``CVLS_NO_ADJ`` -- The function :c:func:`CVodeAdjInit` has not been previously called.
     * ``CVLS_LMEM_NULL`` -- The linear solver has not been initialized with a call to :c:func:`CVodeSetLinearSolverB`.
     * ``CVLS_ILL_INPUT`` -- The parameter ``which`` represented an invalid identifier.


.. c:function:: int CVodeSetLinSysFnBS(void * cvode_mem, int which, CVLsLinSysFnBS linsysBS)

   The function :c:func:`CVodeSetLinSysFnBS` specifies the linear system
   approximation function to be used for the backward problem, in the  case
   where the backward problem depends on the forward sensitivities.

   **Arguments:**
     * ``cvode_mem`` -- pointer to the CVODES memory returned by :c:func:`CVodeCreate`.
     * ``which`` -- represents the identifier of the backward problem.
     * ``linsysBS`` -- user-defined linear system approximation function.

   **Return value:**
     * ``CVLS_SUCCESS`` -- :c:func:`CVodeSetLinSysFnBS` succeeded.
     * ``CVLS_MEM_NULL`` -- ``cvode_mem`` was ``NULL``.
     * ``CVLS_NO_ADJ`` -- The function :c:func:`CVodeAdjInit` has not been previously called.
     * ``CVLS_LMEM_NULL`` -- The linear solver has not been initialized with a call to :c:func:`CVodeSetLinearSolverB`.
     * ``CVLS_ILL_INPUT`` -- The parameter ``which`` represented an invalid identifier.


The function :c:func:`CVodeSetLinearSolutionScalingB` can be used to enable or
disable solution scaling when using a matrix-based linear solver.

.. c:function:: int CVodeSetLinearSolutionScalingB(void * cvode_mem, int which, booleantype onoffB)

   The function :c:func:`CVodeSetLinearSolutionScalingB` enables or disables
   scaling  the linear system solution to account for a change in :math:`\gamma`
   in the linear  system in the backward problem. For more details see
   :numref:`SUNLinsol.CVODES.lagged`.

   **Arguments:**
     * ``cvode_mem`` -- pointer to the CVODES memory block.
     * ``which`` -- represents the identifier of the backward problem.
     * ``onoffB`` -- flag to enable ``SUNTRUE`` or disable ``SUNFALSE`` scaling

   **Return value:**
     * ``CVLS_SUCCESS`` -- The flag value has been successfully set.
     * ``CVLS_MEM_NULL`` -- The ``cvode_mem`` pointer is ``NULL``.
     * ``CVLS_LMEM_NULL`` -- The CVLS linear solver interface has not been initialized.
     * ``CVLS_ILL_INPUT`` -- The attached linear solver is not matrix-based or the linear multistep method type is not BDF.

   **Notes:**
      By default scaling is enabled with matrix-based linear solvers when using
      BDF  methods.


.. c:function:: int CVodeSetJacTimesB(void * cvode_mem, int which, CVLsJacTimesSetupFnB jsetupB, CVLsJacTimesVecFnB jtvB)

   The function :c:func:`CVodeSetJacTimesB` specifies the Jacobian-vector  setup
   and product functions to be used.

   **Arguments:**
     * ``cvode_mem`` -- pointer to the CVODES memory block.
     * ``which`` -- the identifier of the backward problem.
     * ``jtsetupB`` -- user-defined function to set up the Jacobian-vector product. Pass ``NULL`` if no setup is necessary.
     * ``jtvB`` -- user-defined Jacobian-vector product function.

   **Return value:**
     * ``CVLS_SUCCESS`` -- The optional value has been successfully set.
     * ``CVLS_MEM_NULL`` -- ``cvode_mem`` was ``NULL``.
     * ``CVLS_LMEM_NULL`` -- The CVLS linear solver has not been initialized.
     * ``CVLS_NO_ADJ`` -- The function :c:func:`CVodeAdjInit` has not been previously called.
     * ``CVLS_ILL_INPUT`` -- The parameter ``which`` represented an invalid identifier.

   **Notes:**
      The previous routine ``CVSpilsSetJacTimesB`` is now deprecated.


.. c:function:: int CVodeSetJacTimesBS(void * cvode_mem, int which, CVLsJacTimesVecFnBS jtvBS)

   The function :c:func:`CVodeSetJacTimesBS` specifies the Jacobian-vector
   setup and product functions to be used, in the case where the backward
   problem  depends on the forward sensitivities.

   **Arguments:**
     * ``cvode_mem`` -- pointer to the CVODES memory block.
     * ``which`` -- the identifier of the backward problem.
     * ``jtsetupBS`` -- user-defined function to set up the Jacobian-vector product. Pass ``NULL`` if no setup is necessary.
     * ``jtvBS`` -- user-defined Jacobian-vector product function.

   **Return value:**
     * ``CVLS_SUCCESS`` -- The optional value has been successfully set.
     * ``CVLS_MEM_NULL`` -- ``cvode_mem`` was ``NULL``.
     * ``CVLS_LMEM_NULL`` -- The CVLS linear solver has not been initialized.
     * ``CVLS_NO_ADJ`` -- The function :c:func:`CVodeAdjInit` has not been previously called.
     * ``CVLS_ILL_INPUT`` -- The parameter ``which`` represented an invalid identifier.

   **Notes:**
       The previous routine ``CVSpilsSetJacTimesBS`` is now deprecated.


When using the internal difference quotient the user may optionally supply an
alternative right-hand side function for use in the Jacobian-vector product
approximation for the backward problem by calling
:c:func:`CVodeSetJacTimesRhsFnB`. The alternative right-hand side function
should compute a suitable (and differentiable) approximation to the right-hand
side function provided to :c:func:`CVodeInitB` or :c:func:`CVodeInitBS`. For
example, as done in :cite:p:`dorr2010numerical` for a forward integration
without sensitivity analysis, the alternative function may use lagged values
when evaluating a nonlinearity in the right-hand side to avoid differencing a
potentially non-differentiable factor.


.. c:function:: int CVodeSetJacTimesRhsFnB(void * cvode_mem, int which, CVRhsFn jtimesRhsFn)

   The function :c:func:`CVodeSetJacTimesRhsFn` specifies an alternative ODE
   right-hand side function for use in the internal Jacobian-vector product
   difference quotient approximation.

   **Arguments:**
     * ``cvode_mem`` -- pointer to the CVODES memory block.
     * ``which`` -- the identifier of the backward problem.
     * ``jtimesRhsFn`` -- is the CC function which computes the alternative ODE right-hand side function to use in Jacobian-vector product difference quotient approximations.

   **Return value:**
     * ``CVLS_SUCCESS`` -- The optional value has been successfully set.
     * ``CVLS_MEM_NULL`` -- The ``cvode_mem`` pointer is ``NULL``.
     * ``CVLS_LMEM_NULL`` -- The CVLS linear solver has not been initialized.
     * ``CVLS_NO_ADJ`` -- The function :c:func:`CVodeAdjInit` has not been previously called.
     * ``CVLS_ILL_INPUT`` -- The parameter ``which`` represented an invalid identifier or the internal difference quotient approximation is disabled.

   **Notes:**
      The default is to use the right-hand side function provided to
      :c:func:`CVodeInit`  in the internal difference quotient. If the input
      right-hand side function is  ``NULL``, the default is used.  This function
      must be called after the CVLS linear solver interface  has been
      initialized through a call to :c:func:`CVodeSetLinearSolverB`.


.. c:function:: int CVodeSetPreconditionerB(void * cvode_mem, int which, CVLPrecSetupFnB psetupB, CVLsPrecSolveFnB psolveB)

   The function :c:func:`CVodeSetPrecSolveFnB` specifies the preconditioner
   setup and solve functions for the backward integration.

   **Arguments:**
     * ``cvode_mem`` -- pointer to the CVODES memory block.
     * ``which`` -- the identifier of the backward problem.
     * ``psetupB`` -- user-defined preconditioner setup function.
     * ``psolveB`` -- user-defined preconditioner solve function.

   **Return value:**
     * ``CVLS_SUCCESS`` -- The optional value has been successfully set.
     * ``CVLS_MEM_NULL`` -- ``cvode_mem`` was ``NULL``.
     * ``CVLS_LMEM_NULL`` -- The CVLS linear solver has not been initialized.
     * ``CVLS_NO_ADJ`` -- The function :c:func:`CVodeAdjInit` has not been previously called.
     * ``CVLS_ILL_INPUT`` -- The parameter ``which`` represented an invalid identifier.

   **Notes:**
      The ``psetupB`` argument may be ``NULL`` if no setup operation is involved  in the preconditioner.  The previous routine :c:type:`CVSpilsSetPrecSolveFnB` is now deprecated.


.. c:function:: int CVodeSetPreconditionerBS(void * cvode_mem, int which, CVLsPrecSetupFnBS psetupBS, CVLsPrecSolveFnBS psolveBS)

   The function :c:func:`CVodeSetPrecSolveFnBS` specifies the preconditioner
   setup and solve functions for the backward integration, in the case  where
   the backward problem depends on the forward sensitivities.

   **Arguments:**
     * ``cvode_mem`` -- pointer to the CVODES memory block.
     * ``which`` -- the identifier of the backward problem.
     * ``psetupBS`` -- user-defined preconditioner setup function.
     * ``psolveBS`` -- user-defined preconditioner solve function.

   **Return value:**
     * ``CVLS_SUCCESS`` -- The optional value has been successfully set.
     * ``CVLS_MEM_NULL`` -- ``cvode_mem`` was ``NULL``.
     * ``CVLS_LMEM_NULL`` -- The CVLS linear solver has not been initialized.
     * ``CVLS_NO_ADJ`` -- The function :c:func:`CVodeAdjInit` has not been previously called.
     * ``CVLS_ILL_INPUT`` -- The parameter ``which`` represented an invalid identifier.

   **Notes:**
      The ``psetupBS`` argument may be ``NULL`` if no setup operation is
      involved  in the preconditioner.  The previous routine
      :c:type:`CVSpilsSetPrecSolveFnBS` is now deprecated.


.. c:function:: int CVodeSetEpsLinB(void * cvode_mem, int which, realtype eplifacB)

   The function :c:func:`CVodeSetEpsLinB` specifies the factor by  which the
   Krylov linear solver's convergence test constant is reduced  from the
   nonlinear iteration test constant.  This routine can be used in both the
   cases where the backward problem  does and does not depend on the forward
   sensitvities.

   **Arguments:**
     * ``cvode_mem`` -- pointer to the CVODES memory block.
     * ``which`` -- the identifier of the backward problem.
     * ``eplifacB`` -- value of the convergence test constant reduction factor :math:`\geq 0.0`.

   **Return value:**
     * ``CVLS_SUCCESS`` -- The optional value has been successfully set.
     * ``CVLS_MEM_NULL`` -- ``cvode_mem`` was ``NULL``.
     * ``CVLS_LMEM_NULL`` -- The CVLS linear solver has not been initialized.
     * ``CVLS_NO_ADJ`` -- The function :c:func:`CVodeAdjInit` has not been previously called.
     * ``CVLS_ILL_INPUT`` -- The parameter ``which`` represented an invalid identifier, or ``eplifacB`` was negative.

   **Notes:**
      The default value is :math:`0.05`.  Passing a value ``eplifacB = 0.0``
      also indicates using the default value.  The previous routine
      ``CVSpilsSetEpsLinB`` is now deprecated.


.. c:function:: int CVodeSetLSNormFactorB(void * cvode_mem, int which, realtype nrmfac)

   The function :c:func:`CVodeSetLSNormFactor` specifies the factor to use when
   converting from the integrator tolerance (WRMS norm) to the linear solver
   tolerance (L2 norm) for Newton linear system solves e.g.,  ``tol_L2 = fac *
   tol_WRMS``.  This routine can be used in both the cases wherethe backward
   problem  does and does not depend on the forward sensitvities.

   **Arguments:**
     * ``cvode_mem`` -- pointer to the CVODES memory block.
     * ``which`` -- the identifier of the backward problem.
     * ``nrmfac`` -- the norm conversion factor. If ``nrmfac`` is: :math:`> 0` then the provided value is used. :math:`= 0` then the conversion factor is computed using the vector length i.e., ``nrmfac = N_VGetLength(y)`` default. :math:`< 0` then the conversion factor is computed using the vector dot product ``nrmfac = N_VDotProd(v,v)`` where all the entries of ``v`` are one.

   **Return value:**
     * ``CVLS_SUCCESS`` -- The optional value has been successfully set.
     * ``CVLS_MEM_NULL`` -- ``cvode_mem`` was ``NULL``.
     * ``CVLS_LMEM_NULL`` -- The CVLS linear solver has not been initialized.
     * ``CVLS_NO_ADJ`` -- The function :c:func:`CVodeAdjInit` has not been previously called.
     * ``CVLS_ILL_INPUT`` -- The parameter ``which`` represented an invalid identifier.

   **Notes:**
      This function must be called after the CVLS linear solver  interface has
      been initialized through a call to  :c:func:`CVodeSetLinearSolverB`.
      Prior to the introduction of ``N_VGetLength`` in SUNDIALS v5.0.0  (CVODES
      v5.0.0) the value of ``nrmfac`` was computed using the vector  dot product
      i.e., the ``nrmfac < 0`` case.


.. _CVODES.Usage.ADJ.user_callable.optional_output_b:

Optional output functions for the backward problem
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The user of the adjoint module in CVODES has access to any of the optional
output functions described in
:numref:`CVODES.Usage.SIM.optional_output`, both for the main solver
and for the linear solver modules. The first argument of these ``CVodeGet*`` and
``CVode*Get*`` functions is the pointer to the CVODES memory block for the
backward problem. In order to call any of these functions, the user must first
call the following function to obtain this pointer.

.. c:function:: int CVodeGetAdjCVodeBmem(void * cvode_mem, int which)

   The function :c:func:`CVodeGetAdjCVodeBmem` returns a pointer to the CVODES  memory block for the backward problem.

   **Arguments:**
     * ``cvode_mem`` -- pointer to the CVODES memory block created by :c:func:`CVodeCreate`.
     * ``which`` -- the identifier of the backward problem.

   **Return value:**
     * ``void``

   .. warning::
      The user should not modify ``cvode_memB`` in any way.  Optional output calls should pass ``cvode_memB`` as the first argument;  for example, to get the number of integration steps:  ``flag = CVodeGetNumSteps(cvodes_memB, nsteps)``.


To get values of the *forward* solution during a backward integration, use the
following function. The input value of ``t`` would typically be equal to that at
which the backward solution has just been obtained with :c:func:`CVodeGetB`. In
any case, it must be within the last checkpoint interval used by
:c:func:`CVodeB`.

.. c:function:: int CVodeGetAdjY(void * cvode_mem, realtype t, N_Vector y)

   The function :c:func:`CVodeGetAdjY` returns the interpolated value of  the forward solution :math:`y` during a backward integration.

   **Arguments:**
     * ``cvode_mem`` -- pointer to the CVODES memory block created by :c:func:`CVodeCreate`.
     * ``t`` -- value of the independent variable at which :math:`y` is desired input.
     * ``y`` -- forward solution :math:`y(t)`.

   **Return value:**
     * ``CV_SUCCESS`` -- :c:func:`CVodeGetAdjY` was successful.
     * ``CV_MEM_NULL`` -- ``cvode_mem`` was ``NULL``.
     * ``CV_GETY_BADT`` -- The value of ``t`` was outside the current checkpoint interval.

   .. warning::
       The user must allocate space for ``y``.


.. c:function:: int CVodeGetAdjCheckPointsInfo(void * cvode_mem, CVadjCheckPointRec *ckpnt)

   The function :c:func:`CVodeGetAdjCheckPointsInfo` loads an array of ``ncheck+1``  records of type ``CVadjCheckPointRec``.  The user must allocate space for the array ``ckpnt``.

   **Arguments:**
     * ``cvode_mem`` -- pointer to the CVODES memory block created by :c:func:`CVodeCreate`.
     * ``ckpnt`` -- array of ``ncheck+1`` checkpoint records.

   **Return value:**
     * ``void``

   **Notes:**
      The members of each record ``ckpnt[i]`` are:

      * ``ckpnt[i].my_addr`` (``void *``) -- address of current checkpoint in ``cvode_mem->cv_adj_mem``
      * ``ckpnt[i].next_addr`` (``void *``) -- address of next checkpoint
      * ``ckpnt[i].t0`` (``realtype``) -- start of checkpoint interval
      * ``ckpnt[i].t1`` (``realtype``) -- end of checkpoint interval
      * ``ckpnt[i].nstep`` (``long int``) -- step counter at ckeckpoint ``t0``
      * ``ckpnt[i].order`` (``int``) -- method order at checkpoint ``t0``
      * ``ckpnt[i].step`` (``realtype``) -- step size at checkpoint ``t0``


Backward integration of quadrature equations
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Not only the backward problem but also the backward quadrature equations may or
may not depend on the forward sensitivities. Accordingly, either
:c:func:`CVodeQuadInitB` or :c:func:`CVodeQuadInitBS` should be used to allocate
internal memory and to initialize backward quadratures. For any other operation
(extraction, optional input/output, reinitialization, deallocation), the same
function is callable regardless of whether or not the quadratures are
sensitivity-dependent.

.. _CVODES.Usage.ADJ.user_callable.backquad.cvquadinitb:

Backward quadrature initialization functions
""""""""""""""""""""""""""""""""""""""""""""

The function :c:func:`CVodeQuadInitB` initializes and allocates memory for the
backward integration of quadrature equations that do not depend on forward
sensitivities. It has the following form:

.. c:function:: int CVodeQuadInitB(void * cvode_mem, int which, CVQuadRhsFnB rhsQB, N_Vector yQB0)

   The function :c:func:`CVodeQuadInitB` provides required problem specifications,  allocates internal memory, and initializes backward quadrature integration.

   **Arguments:**
     * ``cvode_mem`` -- pointer to the CVODES memory block.
     * ``which`` -- the identifier of the backward problem.
     * ``rhsQB`` -- is the function which computes :math:`fQB`.
     * ``yQB0`` -- is the value of the quadrature variables at ``tB0``.

   **Return value:**
     * ``CV_SUCCESS`` -- The call to :c:func:`CVodeQuadInitB` was successful.
     * ``CV_MEM_NULL`` -- ``cvode_mem`` was ``NULL``.
     * ``CV_NO_ADJ`` -- The function :c:func:`CVodeAdjInit` has not been previously called.
     * ``CV_MEM_FAIL`` -- A memory allocation request has failed.
     * ``CV_ILL_INPUT`` -- The parameter ``which`` is an invalid identifier.


The function :c:func:`CVodeQuadInitBS` initializes and allocates memory for the
backward integration of quadrature equations that depends on the forward
sensitivities.


.. c:function:: int CVodeQuadInitBS(void * cvode_mem, int which, CVQuadRhsFnBS rhsQBS, N_Vector yQBS0)

   The function :c:func:`CVodeQuadInitBS` provides required problem
   specifications,  allocates internal memory, and initializes backward
   quadrature integration.

   **Arguments:**
     * ``cvode_mem`` -- pointer to the CVODES memory block.
     * ``which`` -- the identifier of the backward problem.
     * ``rhsQBS`` -- is the function which computes :math:`fQBS`.
     * ``yQBS0`` -- is the value of the sensitivity-dependent quadrature variables at ``tB0``.

   **Return value:**
     * ``CV_SUCCESS`` -- The call to :c:func:`CVodeQuadInitBS` was successful.
     * ``CV_MEM_NULL`` -- ``cvode_mem`` was ``NULL``.
     * ``CV_NO_ADJ`` -- The function :c:func:`CVodeAdjInit` has not been previously called.
     * ``CV_MEM_FAIL`` -- A memory allocation request has failed.
     * ``CV_ILL_INPUT`` -- The parameter ``which`` is an invalid identifier.


The integration of quadrature equations during the backward phase can be
re-initialized by calling the following function. Before calling
:c:func:`CVodeQuadReInitB` for a new backward problem, call any desired solution
extraction functions ``CVodeGet**`` associated with the previous backward
problem.

.. c:function:: int CVodeQuadReInitB(void * cvode_mem, int which, N_Vector yQB0)

   The function :c:func:`CVodeQuadReInitB` re-initializes the backward quadrature integration.

   **Arguments:**
     * ``cvode_mem`` -- pointer to the CVODES memory block.
     * ``which`` -- the identifier of the backward problem.
     * ``yQB0`` -- is the value of the quadrature variables at ``tB0``.

   **Return value:**
     * ``CV_SUCCESS`` -- The call to :c:func:`CVodeQuadReInitB` was successful.
     * ``CV_MEM_NULL`` -- ``cvode_mem`` was ``NULL``.
     * ``CV_NO_ADJ`` -- The function :c:func:`CVodeAdjInit` has not been previously called.
     * ``CV_MEM_FAIL`` -- A memory allocation request has failed.
     * ``CV_NO_QUAD`` -- Quadrature integration was not activated through a previous call to :c:func:`CVodeQuadInitB`.
     * ``CV_ILL_INPUT`` -- The parameter ``which`` is an invalid identifier.

   **Notes:**
      The function :c:func:`CVodeQuadReInitB` can be called after a call to either  :c:func:`CVodeQuadInitB` or :c:func:`CVodeQuadInitBS`.


.. _CVODES.Usage.ADJ.user_callable.backquad.quad_get_b:

Backward quadrature extraction function
"""""""""""""""""""""""""""""""""""""""

To extract the values of the quadrature variables at the last return time of
:c:func:`CVodeB`, CVODES provides a wrapper for the function
:c:func:`CVodeGetQuad`.


.. c:function:: int CVodeGetQuadB(void * cvode_mem, whichrealtype tret, N_Vector yQB)

   The function :c:func:`CVodeGetQuadB` returns the quadrature solution vector
   after  a successful return from :c:func:`CVodeB`.

   **Arguments:**
     * ``cvode_mem`` -- pointer to the CVODES memory.
     * ``tret`` -- the time reached by the solver output.
     * ``yQB`` -- the computed quadrature vector.

   **Return value:**
     * ``CV_SUCCESS`` -- :c:func:`CVodeGetQuadB` was successful.
     * ``CV_MEM_NULL`` -- ``cvode_mem`` is ``NULL``.
     * ``CV_NO_ADJ`` -- The function :c:func:`CVodeAdjInit` has not been previously called.
     * ``CV_NO_QUAD`` -- Quadrature integration was not initialized.
     * ``CV_BAD_DKY`` -- ``yQB`` was ``NULL``.
     * ``CV_ILL_INPUT`` -- The parameter ``which`` is an invalid identifier.

   .. warning::
      The user must allocate space for ``yQB``.  To obtain the quadratures
      associated with a given backward problem at some  other time within the
      last integration step, first obtain a pointer to the  proper CVODES
      memory structure by calling :c:func:`CVodeGetAdjCVodeBmem`  and then
      use it to call :c:func:`CVodeGetQuadDky`.


.. _CVODES.Usage.ADJ.user_callable.backquad.quad_optional_input_B:

Optional input/output functions for backward quadrature integration
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Optional values controlling the backward integration of quadrature equations can
be changed from their default values through calls to one of the following
functions which are wrappers for the corresponding optional input functions
defined in :numref:`CVODES.Usage.purequad.optional_inputs`. The
user must specify the identifier ``which`` of the backward problem for which the
optional values are specified.

.. code-block:: c

     flag = CVodeSetQuadErrConB(cvode_mem, which, errconQ);
     flag = CVodeQuadSStolerancesB(cvode_mem, which, reltolQ, abstolQ);
     flag = CVodeQuadSVtolerancesB(cvode_mem, which, reltolQ, abstolQ);


Their return value ``flag`` (of type ``int``) can have any of the return values
of its counterparts, but it can also be ``CV_NO_ADJ`` if the function
:c:func:`CVodeAdjInit` has not been previously called or ``CV_ILL_INPUT`` if the
parameter ``which`` was an invalid identifier.

Access to optional outputs related to backward quadrature integration can be
obtained by calling the corresponding ``CVodeGetQuad*`` functions (see
:numref:`CVODES.Usage.purequad.optional_output`). A pointer
``cvode_memB`` to the CVODES memory block for the backward problem, required as
the first argument of these functions, can be obtained through a call to the
functions :c:func:`CVodeGetAdjCVodeBmem`.


.. _CVODES.Usage.ADJ.user_supplied:

User-supplied functions for adjoint sensitivity analysis
--------------------------------------------------------

In addition to the required ODE right-hand side function and any optional
functions for the forward problem, when using the adjoint sensitivity module in
CVODES, the user must supply one function defining the backward problem ODE and,
optionally, functions to supply Jacobian-related information and one or two
functions that define the preconditioner (if an iterative ``SUNLinearSolver``
module is selected) for the backward problem. Type definitions for all these
user-supplied functions are given below.

.. _CVODES.Usage.ADJ.user_supplied.ODErhs_b:

ODE right-hand side for the backward problem
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If the backward problem does not depend on the forward sensitivities, the user
must provide a ``rhsB`` function of type :c:type:`CVRhsFnB` defined as follows:

.. c:type:: int (*CVRhsFnB)(realtype t, N_Vector y, N_Vector yB, N_Vector yBdot, void *user_dataB)

   This function evaluates the right-hand side :math:`f_B(t,y,y_B)` of the
   backward problem  ODE system. This could be either :eq:`CVODES_adj_eqns` or :eq:`CVODES_adj1_eqns`.

   **Arguments:**
     * ``t`` -- is the current value of the independent variable.
     * ``y`` -- is the current value of the forward solution vector.
     * ``yB`` -- is the current value of the backward dependent variable vector.
     * ``yBdot`` -- is the output vector containing the right-hand side :math:`f_B` of the backward ODE problem.
     * ``user_dataB`` -- is a pointer to the same user data passed to :c:func:`CVodeSetUserDataB`.

   **Return value:**
      A :c:type:`CVRhsFnB` should return 0 if successful, a positive value if a recoverable
      error occurred (in which case CVODES will attempt to correct), or a negative
      value if it failed unrecoverably (in which case the integration is halted and
      :c:func:`CVodeB` returns ``CV_RHSFUNC_FAIL``).

   **Notes:**
      Allocation of memory for ``yBdot`` is handled within CVODES.  The ``y``,
      ``yB``, and ``yBdot`` arguments are all  of type ``N_Vector``, but ``yB``
      and ``yBdot`` typically have  different internal representations from
      ``y``. It is the user's  responsibility to access the vector data
      consistently (including the use of the  correct accessor macros from each
      ``N_Vector`` implementation). For the sake of  computational efficiency,
      the vector functions in the two ``N_Vector`` implementations  provided
      with CVODES do not perform any consistency checks with respect to their
      ``N_Vector`` arguments (see :numref:`NVectors`).  The ``user_dataB`` pointer
      is passed to  the user's ``rhsB`` function every time it is called and can
      be the same as the  ``user_data`` pointer used for the forward problem.

   .. warning::
      Before calling the user's ``rhsB`` function, CVODES needs to evaluate
      (through interpolation) the values of the states from the forward
      integration.  If an error occurs in the interpolation, CVODES triggers
      an unrecoverable  failure in the right-hand side function which will
      halt the integration and  :c:func:`CVodeB` will return
      ``CV_RHSFUNC_FAIL``.


.. _CVODES.Usage.ADJ.user_supplied.ODErhs_bs:

ODE right-hand side for the backward problem depending on the forward sensitivities
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If the backward problem does depend on the forward sensitivities, the user must
provide a ``rhsBS`` function of type :c:type:`CVRhsFnBS` defined as follows:

.. c:type:: int (*CVRhsFnBS)(realtype t, N_Vector y, N_Vector *yS, N_Vector yB, N_Vector yBdot, void *user_dataB)

   This function evaluates the right-hand side :math:`f_B(t, y, y_B, s)` of the
   backward problem  ODE system. This could be either :eq:`CVODES_adj_eqns` or
   :eq:`CVODES_adj1_eqns`.

   **Arguments:**
     * ``t`` -- is the current value of the independent variable.
     * ``y`` -- is the current value of the forward solution vector.
     * ``yS`` -- a pointer to an array of ``Ns`` vectors containing the sensitvities of the forward solution.
     * ``yB`` -- is the current value of the backward dependent variable vector.
     * ``yBdot`` -- is the output vector containing the right-hand side.
     * ``user_dataB`` -- is a pointer to user data, same as passed to :c:func:`CVodeSetUserDataB`.

   **Return value:**
      A :c:type:`CVRhsFnBS` should return 0 if successful, a positive value if a
      recoverable error occurred (in which case CVODES will attempt to correct),
      or a negative value if it failed unrecoverably (in which case the
      integration is halted and :c:func:`CVodeB` returns ``CV_RHSFUNC_FAIL``).

   **Notes:**
      Allocation of memory for ``qBdot`` is handled within CVODES.  The ``y``,
      ``yB``, and ``yBdot`` arguments are all of type ``N_Vector``,  but ``yB``
      and ``yBdot`` typically have different internal representations  from
      ``y``.  Likewise for each ``yS[i]``.  It is the user's  responsibility to
      access the vector data consistently (including the use of the  correct
      accessor macros from each ``N_Vector`` implementation). For the sake of
      computational efficiency, the vector functions in the two ``N_Vector``
      implementations  provided with CVODES do not perform any consistency
      checks with respect to their  ``N_Vector`` arguments (see
      :numref:`NVectors`).  The ``user_dataB`` pointer is passed to  the user's
      ``rhsBS`` function every time it is called and can be the same as the
      ``user_data`` pointer used for the forward problem.

   .. warning::
      Before calling the user's ``rhsBS`` function, CVODES needs to evaluate
      (through interpolation) the values of the states from the forward
      integration.  If an error occurs in the interpolation, CVODES triggers
      an unrecoverable  failure in the right-hand side function which will
      halt the integration and  :c:func:`CVodeB` will return
      ``CV_RHSFUNC_FAIL``.


.. _CVODES.Usage.ADJ.user_supplied.rhs_quad_B:

Quadrature right-hand side for the backward problem
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The user must provide an ``fQB`` function of type :c:type:`CVQuadRhsFnB` defined
by

.. c:type:: int (*CVQuadRhsFnB)(realtype t, N_Vector y, N_Vector yB, N_Vector qBdot, void *user_dataB)

   This function computes the quadrature equation right-hand side for the
   backward problem.

   **Arguments:**
     * ``t`` -- is the current value of the independent variable.
     * ``y`` -- is the current value of the forward solution vector.
     * ``yB`` -- is the current value of the backward dependent variable vector.
     * ``qBdot`` -- is the output vector containing the right-hand side ``fQB`` of the backward quadrature equations.
     * ``user_dataB`` -- is a pointer to user data, same as passed to :c:func:`CVodeSetUserDataB`.

   **Return value:**
      A :c:type:`CVQuadRhsFnB` should return 0 if successful, a positive value
      if a recoverable error occurred (in which case CVODES will attempt to
      correct), or a negative value if it failed unrecoverably (in which case
      the integration is halted and :c:func:`CVodeB` returns
      ``CV_QRHSFUNC_FAIL``).

   **Notes:**
      Allocation of memory for ``rhsvalBQ`` is handled within CVODES.  The
      ``y``, ``yB``, and ``qBdot`` arguments are all of type ``N_Vector``,  but
      they typically do not all have the same representation. It is the user's
      responsibility to access the vector data consistently (including the use
      of the  correct accessor macros from each ``N_Vector`` implementation).
      For the sake of  computational efficiency, the vector functions in the two
      ``N_Vector`` implementations  provided with CVODES do not perform any
      consistency checks with repsect to their  ``N_Vector`` arguments
      (see :numref:`NVectors`).  The ``user_dataB`` pointer is passed to the user's
      ``fQB`` function every time  it is called and can be the same as the
      ``user_data`` pointer used for the forward problem.

   .. warning::
      Before calling the user's ``fQB`` function, CVODES needs to evaluate
      (through interpolation) the values of the states from the forward
      integration.  If an error occurs in the interpolation, CVODES triggers
      an unrecoverable  failure in the quadrature right-hand side function
      which will halt the integration and  :c:func:`CVodeB` will return
      ``CV_QRHSFUNC_FAIL``.


.. _CVODES.Usage.ADJ.user_supplied.rhs_quad_sens_B:

Sensitivity-dependent quadrature right-hand side for the backward problem
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The user must provide an ``fQBS`` function of type :c:type:`CVQuadRhsFnBS`
defined by

.. c:type:: int (*CVQuadRhsFnBS)(realtype t, N_Vector y, N_Vector *yS, N_Vector yB, N_Vector qBdot, void *user_dataB)

   This function computes the quadrature equation right-hand side for the
   backward problem.

   **Arguments:**
     * ``t`` -- is the current value of the independent variable.
     * ``y`` -- is the current value of the forward solution vector.
     * ``yS`` -- a pointer to an array of ``Ns`` vectors continaing the sensitvities of the forward solution.
     * ``yB`` -- is the current value of the backward dependent variable vector.
     * ``qBdot`` -- is the output vector containing the right-hand side ``fQBS`` of the backward quadrature equations.
     * ``user_dataB`` -- is a pointer to user data, same as passed to :c:func:`CVodeSetUserDataB`.

   **Return value:**
      A :c:type:`CVQuadRhsFnBS` should return 0 if successful, a positive value if a recoverable
      error occurred (in which case CVODES will attempt to correct), or a negative
      value if it failed unrecoverably (in which case the integration is halted and
      :c:func:`CVodeB` returns ``CV_QRHSFUNC_FAIL``).

   **Notes:**
      Allocation of memory for ``qBdot`` is handled within CVODES.  The ``y``,
      ``yS``, and ``qBdot`` arguments are all of type ``N_Vector``,  but they
      typically do not all have the same internal representation.  Likewise for
      each ``yS[i]``.  It is the user's  responsibility to access the vector
      data consistently (including the use of the  correct accessor macros from
      each ``N_Vector`` implementation). For the sake of  computational
      efficiency, the vector functions in the two ``N_Vector`` implementations
      provided with CVODES do not perform any consistency checks with repsect to
      their  ``N_Vector`` arguments (see :numref:`NVectors`).  The
      ``user_dataB`` pointer is passed to the user's ``fQBS`` function every
      time  it is called and can be the same as the ``user_data`` pointer used
      for the forward problem.

   .. warning::
      Before calling the user's ``fQBS`` function, CVODES needs to evaluate
      (through interpolation) the values of the states from the forward
      integration.  If an error occurs in the interpolation, CVODES triggers
      an unrecoverable  failure in the quadrature right-hand side function
      which will halt the integration and  :c:func:`CVodeB` will return
      ``CV_QRHSFUNC_FAIL``.


.. _CVODES.Usage.ADJ.user_supplied.jacFn_b:

Jacobian construction for the backward problem (matrix-based linear solvers)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If a matrix-based linear solver module is used for the backward problem (i.e., a
non-``NULL`` ``SUNMatrix`` object was supplied to
:c:func:`CVodeSetLinearSolverB`), the user may provide a function of type
:c:type:`CVLsJacFnB` or :c:type:`CVLsJacFnBS`, defined as follows:

.. c:type:: int (*CVLsJacFnB)(realtype t, N_Vector y, N_Vector yB, N_Vector fyB, SUNMatrix JacB, void *user_dataB, N_Vector tmp1B, N_Vector tmp2B, N_Vector tmp3B)

   This function computes the Jacobian of the backward problem (or an
   approximation  to it).

   **Arguments:**
     * ``t`` -- is the current value of the independent variable.
     * ``y`` -- is the current value of the forward solution vector.
     * ``yB`` -- is the current value of the backward dependent variable vector.
     * ``fyB`` -- is the current value of the backward right-hand side function :math:`f_B`.
     * ``JacB`` -- is the output approximate Jacobian matrix.
     * ``user_dataB`` -- is a pointer to the same user data passed to :c:func:`CVodeSetUserDataB`.
     * ``tmp1B``, ``tmp2B``, ``tmp3B`` -- are pointers to memory allocated for variables of type ``N_Vector`` which can be used by the :c:type:`CVLsJacFnB` function as temporary storage or work space.

   **Return value:**
      A :c:type:`CVLsJacFnB` should return 0 if successful, a positive value if a recoverable
      error occurred (in which case CVODES will attempt to correct, while CVLS sets
      ``last_flag`` to ``CVLS_JACFUNC_RECVR``), or a negative
      value if it failed unrecoverably (in which case the integration is halted, :c:func:`CVodeB`
      returns ``CV_LSETUP_FAIL`` and CVLS sets ``last_flag`` to
      ``CVLS_JACFUNC_UNRECVR``).

   **Notes:**
      A user-supplied Jacobian function must load the  matrix ``JacB`` with an
      approximation to the Jacobian matrix  at the point ``(t, y, yB)``,
      where ``y`` is the solution  of the original IVP at time ``tt``, and
      ``yB`` is the solution of the  backward problem at the same time.
      Information regarding the structure of the specific ``SUNMatrix``
      structure (e.g. number of rows, upper/lower bandwidth, sparsity  type) may
      be obtained through using the implementation-specific  ``SUNMatrix``
      interface functions (see :numref:`SUNMatrix` for details).  With direct linear
      solvers (i.e., linear solvers with type  ``SUNLINEARSOLVER_DIRECT``), the
      Jacobian matrix :math:`J(t,y)` is zeroed out  prior to calling the
      user-supplied Jacobian function so only nonzero elements  need to be
      loaded into ``JacB``.

   .. warning::
      Before calling the user's :c:type:`CVLsJacFnB`, CVODES needs to
      evaluate  (through interpolation) the values of the states from the
      forward integration.  If an error occurs in the interpolation, CVODES
      triggers an unrecoverable  failure in the Jacobian function which will
      halt the integration  (:c:func:`CVodeB` returns ``CV_LSETUP_FAIL`` and
      CVLS sets ``last_flag`` to  ``CVLS_JACFUNC_UNRECVR``).  The previous
      function type :c:type:`CVDlsJacFnB` is identical to
      :c:type:`CVLsJacFnB`, and may still be used for backward-compatibility.
      However, this will be deprecated in future releases, so we recommend
      that users transition to the new function type name soon.


.. c:type:: int (*CVLsJacFnBS)(realtype t, N_Vector y, N_Vector *yS, N_Vector yB, N_Vector fyB, SUNMatrix JacB, void *user_dataB, N_Vector tmp1B, N_Vector tmp2B, N_Vector tmp3B)

   This function computes the Jacobian of the backward problem (or an
   approximation to it), in the case where the backward problem depends on the
   forward sensitivities.

   **Arguments:**
     * ``t`` -- is the current value of the independent variable.
     * ``y`` -- is the current value of the forward solution vector.
     * ``yS`` -- a pointer to an array of ``Ns`` vectors containing the sensitvities of the forward solution.
     * ``yB`` -- is the current value of the backward dependent variable vector.
     * ``fyB`` -- is the current value of the backward right-hand side function :math:`f_B`.
     * ``JacB`` -- is the output approximate Jacobian matrix.
     * ``user_dataB`` -- is a pointer to the same user data passed to :c:func:`CVodeSetUserDataB`.
     * ``tmp1B``, ``tmp2B``, ``tmp3B`` -- are pointers to memory allocated for variables of type ``N_Vector`` which can be used by the :c:type:`CVLsLinSysFnBS` function as temporary storage or work space.

   **Return value:**
      A :c:type:`CVLsJacFnBS` should return 0 if successful, a positive value if a recoverable
      error occurred (in which case CVODES will attempt to correct, while CVLS sets
      ``last_flag`` to ``CVLS_JACFUNC_RECVR``), or a negative
      value if it failed unrecoverably (in which case the integration is halted, :c:func:`CVodeB`
      returns ``CV_LSETUP_FAIL`` and CVLS sets ``last_flag`` to
      ``CVLS_JACFUNC_UNRECVR``).

   **Notes:**
      A user-supplied Jacobian function must load the  matrix ``JacB`` with an
      approximation to the Jacobian matrix at the point
      ``(t, y, yS, yB)``, where ``y`` is the solution of the original
      IVP at time ``tt``, ``yS`` is the vector of forward sensitivities at time
      ``tt``,  and ``yB`` is the solution of the backward problem at the same
      time.  Information regarding the structure of the specific ``SUNMatrix``
      structure (e.g. number of rows, upper/lower bandwidth, sparsity  type) may
      be obtained through using the implementation-specific  ``SUNMatrix``
      interface functions (see :numref:`SUNMatrix`).  With direct linear solvers
      (i.e., linear solvers with type  ``SUNLINEARSOLVER_DIRECT``, the Jacobian
      matrix :math:`J(t,y)` is zeroed out prior  to calling the user-supplied
      Jacobian function so only nonzero elements need  to be loaded into
      ``JacB``.

   .. warning::
      Before calling the user's :c:type:`CVLsJacFnBS`, CVODES needs to
      evaluate  (through interpolation) the values of the states from the
      forward integration.  If an error occurs in the interpolation, CVODES
      triggers an unrecoverable  failure in the Jacobian function which will
      halt the integration  (:c:func:`CVodeB` returns ``CV_LSETUP_FAIL`` and
      CVLS sets ``last_flag`` to  ``CVLS_JACFUNC_UNRECVR``).  The previous
      function type :c:type:`CVDlsJacFnBS` is identical to
      :c:type:`CVLsJacFnBS`, and may still be used for
      backward-compatibility.  However, this will be deprecated in future
      releases, so we recommend  that users transition to the new function
      type name soon.



.. _CVODES.Usage.ADJ.user_supplied.linsysFn_b:

Linear system construction for the backward problem (matrix-based linear solvers)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

With matrix-based linear solver modules, as an alternative to optionally
supplying a function for evaluating the Jacobian of the ODE right-hand side
function, the user may optionally supply a function of type :c:type:`CVLsLinSysFnB` or
:c:type:`CVLsLinSysFnBS` for evaluating the linear system, :math:`M_B = I - \gamma_B
J_B` (or an approximation of it) for the backward problem.

.. c:type:: int (*CVLsLinSysFnB)(realtype t, N_Vector y, N_Vector yB, N_Vector fyB, SUNMatrix AB, booleantype jokB, booleantype *jcurB, realtype gammaB, void *user_dataB, N_Vector tmp1B, N_Vector tmp2B, N_Vector tmp3B);

   This function computes the linear system of the backward problem (or an
   approximation to it).

   **Arguments:**
     * ``t`` -- is the current value of the independent variable.
     * ``y`` -- is the current value of the forward solution vector.
     * ``yB`` -- is the current value of the backward dependent variable vector.
     * ``fyB`` -- is the current value of the backward right-hand side function :math:`f_B`.
     * ``AB`` -- is the output approximate linear system matrix.
     * ``jokB`` -- is an input flag indicating whether Jacobian-related data needs to be recomputed (``jokB = SUNFALSE``) or informtion saved from a previous information can be safely used (``jokB = SUNTRUE``).
     * ``jcurB`` -- is an output flag which must be set to ``SUNTRUE`` if Jacobian-related data was recomputed or ``SUNFALSE`` otherwise.
     * ``gammaB`` -- is the scalar appearing in the matrix :math:`M_B = I - \gamma_B J_B`.
     * ``user_dataB`` -- is a pointer to the same user data passed to :c:func:`CVodeSetUserDataB`.
     * ``tmp1B``, ``tmp2B``, ``tmp3B`` -- are pointers to memory allocated for variables of type ``N_Vector`` which can be used by the :c:type:`CVLsLinSysFnB` function as temporary storage or work space.

   **Return value:**
      A :c:type:`CVLsLinSysFnB` should return ``0`` if successful, a positive value if a
      recoverable error occurred (in which case CVODES will attempt to correct,
      while CVLS sets ``last_flag`` to ``CVLS_JACFUNC_RECVR``), or a
      negative value if it failed unrecoverably (in which case the integration is
      halted, :c:func:`CVodeB` returns ``CV_LSETUP_FAIL`` and CVLS sets
      ``last_flag`` to ``CVLS_JACFUNC_UNRECVR``).

   **Notes:**
      A user-supplied linear system function must load the matrix ``AB`` with an
      approximation to the linear system matrix at the point ``(t, y, yB)``,
      where ``y`` is the solution of the original IVP at time ``tt``,
      and ``yB`` is the solution of the backward problem at the same time.

   .. warning::
      Before calling the user's :c:type:`CVLsLinSysFnB`, CVODES needs
      to  evaluate (through interpolation) the values of the states from the
      forward  integration. If an error occurs in the interpolation,
      CVODES triggers an  unrecoverable failure in the linear system
      function which will halt the  integration (:c:func:`CVodeB` returns
      ``CV_LSETUP_FAIL`` and CVLS sets  ``last_flag`` to
      ``CVLS_JACFUNC_UNRECVR``).


.. c:type:: int (*CVLsLinSysFnBS)(realtype t, N_Vector y, N_Vector* yS, N_Vector yB, N_Vector fyB, SUNMatrix AB, booleantype jokB, booleantype *jcurB, realtype gammaB, void *user_dataB, N_Vector tmp1B, N_Vector tmp2B, N_Vector tmp3B);

   This function computes the linear system of the backward problem (or an
   approximation to it), in the case where the backward problem depends on the
   forward sensitivities.

   **Arguments:**
     * ``t`` -- is the current value of the independent variable.
     * ``y`` -- is the current value of the forward solution vector.
     * ``yS`` -- a pointer to an array of ``Ns`` vectors containing the sensitivities of the forward solution.
     * ``yB`` -- is the current value of the backward dependent variable vector.
     * ``fyB`` -- is the current value of the backward right-hand side function :math:`f_B`.
     * ``AB`` -- is the output approximate linear system matrix.
     * ``jokB`` -- is an input flag indicating whether Jacobian-related data needs to be recomputed (``jokB = SUNFALSE``) or informtion saved from a previous information can be safely used (``jokB = SUNTRUE``).
     * ``jcurB`` -- is an output flag which must be set to ``SUNTRUE`` if Jacobian-related data was recomputed or ``SUNFALSE`` otherwise.
     * ``gammaB`` -- is the scalar appearing in the matrix
     * ``user_dataB`` -- is a pointer to the same user data passed to :c:func:`CVodeSetUserDataB`.
     * ``tmp1B``, ``tmp2B``, ``tmp3B`` -- are pointers to memory allocated for variables of type ``N_Vector`` which can be used by the :c:type:`CVLsLinSysFnBS` function as temporary storage or work space.

   **Return value:**
      A :c:type:`CVLsLinSysFnBS` should return ``0`` if successful, a positive
      value if a recoverable error occurred (in which case CVODES will attempt
      to correct, while CVLS sets ``last_flag`` to ``CVLS_JACFUNC_RECVR``), or a
      negative value if it failed unrecoverably (in which case the integration
      is halted, :c:func:`CVodeB` returns ``CV_LSETUP_FAIL`` and CVLS sets
      ``last_flag`` to ``CVLS_JACFUNC_UNRECVR``).

   **Notes:**
      A user-supplied linear system function must load the matrix ``AB`` with an
      approximation to the linear system matrix at the point ``(t, y, yS, yB)``,
      where ``y`` is the solution of the original IVP at time
      ``tt``, ``yS`` is the vector of forward sensitivities at time ``t``, and
      ``yB`` is the solution of the backward problem at the same time.

   .. warning::
      Before calling the user's :c:type:`CVLsLinSysFnBS`, CVODES needs to
      evaluate (through interpolation) the values of the states from the
      forward  integration.  If an error occurs in the interpolation, CVODES
      triggers an  unrecoverable failure in the linear system function which
      will halt the  integration (:c:func:`CVodeB` returns ``CV_LSETUP_FAIL``
      and CVLS sets  ``last_flag`` to ``CVLS_JACFUNC_UNRECVR``).



.. _CVODES.Usage.ADJ.user_supplied.jtimesv_b:

Jacobian-vector product for the backward problem (matrix-free linear solvers)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If a matrix-free linear solver is to be used for the backward problem (i.e., a
``NULL``-valued ``SUNMatrix`` was supplied to :c:func:`CVodeSetLinearSolverB` in
the steps described in :numref:`CVODES.Usage.ADJ.skeleton_sim`, the user may
provide a function of type :c:type:`CVLsJacTimesVecFnB` or
:c:type:`CVLsJacTimesVecFnBS` in the following form, to compute matrix-vector
products :math:`Jv`. If such a function is not supplied, the default is a
difference quotient approximation to these products.

.. c:type:: int (*CVLsJacTimesVecFnB)(N_Vector vB, N_Vector JvB, realtype t, N_Vector y, N_Vector yB, N_Vector fyB, void *jac_dataB, N_Vector tmpB);

   This function computes the action of the Jacobian ``JB`` for  the backward
   problem on a given vector ``vB``.

   **Arguments:**
     * ``vB`` -- is the vector by which the Jacobian must be multiplied to the right.
     * ``JvB`` -- is the computed output vector ``JB*vB``.
     * ``t`` -- is the current value of the independent variable.
     * ``y`` -- is the current value of the forward solution vector.
     * ``yB`` -- is the current value of the backward dependent variable vector.
     * ``fyB`` -- is the current value of the backward right-hand side function :math:`f_B`.
     * ``user_dataB`` -- is a pointer to the same user data passed to :c:func:`CVodeSetUserDataB`.
     * ``tmpB`` -- is a pointer to memory allocated for a variable of type ``N_Vector`` which can be used by :c:type:`CVLsJacTimesVecFnB` as temporary storage or work space.

   **Return value:**
      The return value of a function of type :c:type:`CVLsJacTimesVecFnB` should be
      if successful or nonzero if an error was encountered, in which case
      the integration is halted.

   **Notes:**
      A user-supplied Jacobian-vector product function must load the  vector
      ``JvB`` with the product of the Jacobian of the backward  problem at the
      point ``(t, y, yB)`` and the vector ``vB``.  Here, ``y`` is the
      solution of the original IVP at time ``t`` and  ``yB`` is the solution of
      the backward problem at the same time.  The rest of the arguments are
      equivalent to those passed to a function of type
      :c:type:`CVLsJacTimesVecFn`.  If the backward problem is the adjoint of
      :math:`{\dot y} = f(t, y)`, then this  function is to compute
      :math:`-({\partial f}/{\partial y_i})^T v_B`.  The previous function type
      :c:type:`CVSpilsJacTimesVecFnB` is deprecated.


.. c:type:: int (*CVLsJacTimesVecFnBS)(N_Vector vB, N_Vector JvB, realtype t, N_Vector y, N_Vector *yS, N_Vector yB, N_Vector fyB, void *user_dataB, N_Vector tmpB);

   This function computes the action of the Jacobian ``JB`` for  the backward
   problem on a given vector ``vB``, in the case where  the backward problem
   depends on the forward sensitivities.

   **Arguments:**
     * ``vB`` -- is the vector by which the Jacobian must be multiplied to the right.
     * ``JvB`` -- is the computed output vector ``JB*vB``.
     * ``t`` -- is the current value of the independent variable.
     * ``y`` -- is the current value of the forward solution vector.
     * ``yS`` -- is a pointer to an array containing the forward sensitivity vectors.
     * ``yB`` -- is the current value of the backward dependent variable vector.
     * ``fyB`` -- is the current value of the backward right-hand side function :math:`f_B`.
     * ``user_dataB`` -- is a pointer to the same user data passed to :c:func:`CVodeSetUserDataB`.
     * ``tmpB`` -- is a pointer to memory allocated for a variable of type ``N_Vector`` which can be used by :c:type:`CVLsJacTimesVecFnB` as temporary storage or work space.

   **Return value:**
      The return value of a function of type :c:type:`CVLsJacTimesVecFnBS` should be
      if successful or nonzero if an error was encountered, in which case
      the integration is halted.

   **Notes:**
      A user-supplied Jacobian-vector product function must load the vector
      ``JvB``  with the product of the Jacobian of the backward problem  at the
      point ``(t, y, yB)`` and the vector ``vB``.  Here, ``y`` is the
      solution of the original IVP at time ``t`` and  ``yB`` is the solution of
      the backward problem at the same time.  The rest of the arguments are
      equivalent to those passed to a function of type
      :c:type:`CVLsJacTimesVecFn`.  The previous function type
      :c:type:`CVSpilsJacTimesVecFnBS` is deprecated.


.. _CVODES.Usage.ADJ.user_supplied.jactimesvecsetup_b:

Jacobian-vector product setup for the backward problem (matrix-free linear solvers)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If the user’s Jacobian-times-vector routine requires that any Jacobian-related
data be preprocessed or evaluated, then this needs to be done in a user-supplied
function of type :c:type:`CVLsJacTimesSetupFnB` or
:c:type:`CVLsJacTimesSetupFnBS`, defined as follows:

.. c:type:: int (*CVLsJacTimesSetupFnB)(realtype t, N_Vector y, N_Vector yB, N_Vector fyB, void *user_dataB)

   This function preprocesses and/or evaluates Jacobian data needed  by the
   Jacobian-times-vector routine for the backward problem.

   **Arguments:**
     * ``t`` -- is the current value of the independent variable.
     * ``y`` -- is the current value of the dependent variable vector, :math:`y(t)`.
     * ``yB`` -- is the current value of the backward dependent variable vector.
     * ``fyB`` -- is the current value of the right-hand-side for the backward problem.
     * ``user_dataB`` -- is a pointer to user data :c:func:`CVodeSetUserDataB`.

   **Return value:**
      The value returned by the Jacobian-vector setup function
      should be  if successful, positive for a recoverable error (in
      which case the step will be retried), or negative for an
      unrecoverable error (in which case the integration is halted).

   **Notes:**
      Each call to the Jacobian-vector setup function is preceded by a call to
      the backward problem residual user function with the same  ``(t,y, yB)``
      arguments.  Thus, the setup function can use any auxiliary data that is
      computed  and saved during the evaluation of the right-hand-side function.
      If the user's :c:type:`CVLsJacTimesVecFnB` function uses difference
      quotient approximations, it may need to access quantities not in the call
      list. These include the current stepsize, the error weights, etc.  To
      obtain these, the user will need to add a pointer to ``cvode_mem``  to
      ``user_dataB`` and then use the ``CVGet*`` functions described in
      :numref:`CVODES.Usage.SIM.optional_output`. The unit
      roundoff can be accessed as  ``UNIT_ROUNDOFF`` defined in
      ``sundials_types.h``.  The previous function type
      :c:type:`CVSpilsJacTimesSetupFnB` is identical  to
      :c:type:`CVLsJacTimesSetupFnB`, and may still be used for
      backward-compatibility.  However, this will be deprecated in future
      releases, so we recommend that users transition to the new function  type
      name soon.


.. c:type:: int (*CVLsJacTimesSetupFnBS)(realtype t, N_Vector y, N_Vector *yS, N_Vector yB, N_Vector fyB, void *user_dataB)

   This function preprocesses and/or evaluates Jacobian data needed  by the
   Jacobian-times-vector routine for the backward problem, in the case that  the
   backward problem depends on the forward sensitivities.

   **Arguments:**
     * ``t`` -- is the current value of the independent variable.
     * ``y`` -- is the current value of the dependent variable vector, :math:`y(t)`.
     * ``yS`` -- a pointer to an array of ``Ns`` vectors containing the sensitvities of the forward solution.
     * ``yB`` -- is the current value of the backward dependent variable vector.
     * ``fyB`` -- is the current value of the right-hand-side function for the backward problem.
     * ``user_dataB`` -- is a pointer to the same user data provided to :c:func:`CVodeSetUserDataB`.

   **Return value:**
      The value returned by the Jacobian-vector setup function
      should be  if successful, positive for a recoverable error (in
      which case the step will be retried), or negative for an
      unrecoverable error (in which case the integration is halted).

   **Notes:**
      Each call to the Jacobian-vector setup function is preceded by a call to
      the backward problem residual user function with the same  ``(t,y, yS,
      yB)`` arguments.  Thus, the setup function can use any auxiliary data that
      is computed  and saved during the evaluation of the right-hand-side
      function.  If the user's :c:type:`CVLsJacTimesVecFnBS` function uses
      difference quotient  approximations, it may need to access quantities not
      in the call list. These include the current stepsize, the error weights,
      etc.  To obtain these, the user will need to add a pointer to
      ``cvode_mem``  to ``user_dataB`` and then use the ``CVGet*`` functions
      described in  :numref:`CVODES.Usage.SIM.optional_output`.
      The unit roundoff can be accessed as  ``UNIT_ROUNDOFF`` defined in
      ``sundials_types.h``.  The previous function type
      :c:type:`CVSpilsJacTimesSetupFnBS` is identical  to
      :c:type:`CVLsJacTimesSetupFnBS`, and may still be used for
      backward-compatibility.  However, this will be deprecated in future
      releases, so we recommend that users transition to the new function  type
      name soon.


.. _CVODES.Usage.ADJ.user_supplied.psolve_b:

Preconditioner solve for the backward problem (iterative linear solvers)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If a user-supplied preconditioner is to be used with a ``SUNLinearSolver``
solver module, then the user must provide a function to solve the linear system
:math:`Pz = r`, where :math:`P` may be either a left or a right preconditioner
matrix. Here :math:`P` should approximate (at least crudely) the matrix
:math:`M_B = I - \gamma_B J_B`, where :math:`J_B = \partial f_B/ \partial y_B`.
If preconditioning is done on both sides, the product of the two preconditioner
matrices should approximate :math:`M_B`. This function must be of one of the
following two types:

.. c:type:: int (*CVLsPrecSolveFnB)(realtype t, N_Vector y, N_Vector yB, N_Vector fyB, N_Vector rvecB, N_Vector zvecB, realtype gammaB, realtype deltaB, void *user_dataB)

   This function solves the preconditioning system :math:`Pz = r` for the
   backward problem.

   **Arguments:**
     * ``t`` -- is the current value of the independent variable.
     * ``y`` -- is the current value of the forward solution vector.
     * ``yB`` -- is the current value of the backward dependent variable vector.
     * ``fyB`` -- is the current value of the backward right-hand side function :math:`f_B`.
     * ``rvecB`` -- is the right-hand side vector ``r`` of the linear system to be solved.
     * ``zvecB`` -- is the computed output vector.
     * ``gammaB`` -- is the scalar appearing in the matrix, :math:`M_B = I - \gamma_B J_B`.
     * ``deltaB`` -- is an input tolerance to be used if an iterative method is employed in the solution.
     * ``user_dataB`` -- is a pointer to the same user data passed to :c:func:`CVodeSetUserDataB`.

   **Return value:**
      The return value of a preconditioner solve function for the backward
      problem should be  if successful, positive for a recoverable error (in
      which case the step will be retried), or negative for an unrecoverable
      error (in which case the integration is halted).

   **Notes:**
      The previous function type :c:type:`CVSpilsPrecSolveFnB` is deprecated.


.. c:type:: int (*CVLsPrecSolveFnBS)(realtype t, N_Vector y, N_Vector *yS, N_Vector yB, N_Vector fyB, N_Vector rvecB, N_Vector zvecB, realtype gammaB, realtype deltaB, void *user_dataB)

   This function solves the preconditioning system :math:`Pz = r` for the backward problem,  in the case where the backward problem depends on the forward sensitivities.

   **Arguments:**
     * ``t`` -- is the current value of the independent variable.
     * ``y`` -- is the current value of the forward solution vector.
     * ``yS`` -- is a pointer to an array containing the forward sensitivity vectors.
     * ``yB`` -- is the current value of the backward dependent variable vector.
     * ``fyB`` -- is the current value of the backward right-hand side function :math:`f_B`.
     * ``rvecB`` -- is the right-hand side vector ``r`` of the linear system to be solved.
     * ``zvecB`` -- is the computed output vector.
     * ``gammaB`` -- is the scalar appearing in the matrix, :math:`M_B = I - \gamma_B J_B`.
     * ``deltaB`` -- is an input tolerance to be used if an iterative method is employed in the solution.
     * ``user_dataB`` -- is a pointer to the same user data passed to :c:func:`CVodeSetUserDataB`.

   **Return value:**
      The return value of a preconditioner solve function for the backward
      problem should be  if successful,
      positive for a recoverable error (in which case the step will be retried), or
      negative for an unrecoverable error (in which case the integration is halted).


   **Notes:**
      The previous function type :c:type:`CVSpilsPrecSolveFnBS` is identical to  :c:type:`CVLsPrecSolveFnBS`, and may still be used for backward-compatibility.  However, this will be deprecated in future releases, so we recommend  that users transition to the new function type name soon.


.. _CVODES.Usage.ADJ.user_supplied.psetup_b:

Preconditioner setup for the backward problem (iterative linear solvers)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If the user’s preconditioner requires that any Jacobian-related data be
preprocessed or evaluated, then this needs to be done in a user-supplied
function of one of the following two types:

.. c:type:: int (*CVLsPrecSetupFnB)(realtype t, N_Vector y, N_Vector yB, N_Vector fyB, booleantype jokB, booleantype *jcurPtrB, realtype gammaB, void *user_dataB)

   This function preprocesses and/or evaluates Jacobian-related data needed  by
   the preconditioner for the backward problem.

   **Arguments:**
     * ``t`` -- is the current value of the independent variable.
     * ``y`` -- is the current value of the forward solution vector.
     * ``yB`` -- is the current value of the backward dependent variable vector.
     * ``fyB`` -- is the current value of the backward right-hand side function :math:`f_B`.
     * ``jokB`` -- is an input flag indicating whether Jacobian-related data needs to be recomputed (``jokB = SUNFALSE``) or information saved from a previous invokation can be safely used (``jokB = SUNTRUE``).
     * ``jcurPtr`` -- is an output flag which must be set to ``SUNTRUE`` if Jacobian-related data was recomputed or ``SUNFALSE`` otherwise.
     * ``gammaB`` -- is the scalar appearing in the matrix :math:`M_B = I - \gamma_B J_B`.
     * ``user_dataB`` -- is a pointer to the same user data passed to :c:func:`CVodeSetUserDataB`.

   **Return value:**
      The return value of a preconditioner setup function for the backward
      problem should be  if successful, positive for a recoverable error (in
      which case the step will be retried), or negative for an unrecoverable
      error (in which case the integration is halted).

   **Notes:**
      The previous function type :c:type:`CVSpilsPrecSetupFnB` is identical to
      :c:type:`CVLsPrecSetupFnB`, and may still be used for
      backward-compatibility.  However, this will be deprecated in future
      releases, so we recommend  that users transition to the new function type
      name soon.


.. c:type:: int (*CVLsPrecSetupFnBS)(realtype t, N_Vector y, N_Vector *yS, N_Vector yB, N_Vector fyB, booleantype jokB, booleantype *jcurPtrB, realtype gammaB, void *user_dataB)

   This function preprocesses and/or evaluates Jacobian-related data needed  by
   the preconditioner for the backward problem, in the case where the  backward
   problem depends on the forward sensitivities.

   **Arguments:**
     * ``t`` -- is the current value of the independent variable.
     * ``y`` -- is the current value of the forward solution vector.
     * ``yS`` -- is a pointer to an array containing the forward sensitivity vectors.
     * ``yB`` -- is the current value of the backward dependent variable vector.
     * ``fyB`` -- is the current value of the backward right-hand side function :math:`f_B`.
     * ``jokB`` -- is an input flag indicating whether Jacobian-related data needs to be recomputed (``jokB = SUNFALSE``) or information saved from a previous invokation can be safely used (``jokB = SUNTRUE``).
     * ``jcurPtr`` -- is an output flag which must be set to ``SUNTRUE`` if Jacobian-related data was recomputed or ``SUNFALSE`` otherwise.
     * ``gammaB`` -- is the scalar appearing in the matrix :math:`M_B = I - \gamma_B J_B`.
     * ``user_dataB`` -- is a pointer to the same user data passed to :c:func:`CVodeSetUserDataB`.

   **Return value:**
      The return value of a preconditioner setup function for the backward
      problem should be  if successful, positive for a recoverable error (in
      which case the step will be retried), or negative for an unrecoverable
      error (in which case the integration is halted).

   **Notes:**
      The previous function type :c:type:`CVSpilsPrecSetupFnBS` is deprecated.


Using CVODES preconditioner modules for the backward problem
------------------------------------------------------------

As on the forward integration phase, the efficiency of Krylov iterative methods
for the solution of linear systems can be greatly enhanced through
preconditioning. Both preconditioner modules provided with SUNDIALS, the serial
banded preconditioner CVBANDPRE and the parallel band-block-diagonal
preconditioner module CVBBDPRE, provide interface functions through which they
can be used on the backward integration phase.

Using the banded preconditioner CVBANDPRE
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The adjoint module in CVODES offers an interface to the banded preconditioner
module CVBANDPRE described in section
:numref:`CVODES.Usage.SIM.precond.cvbandpre`. This preconditioner, usable only in
a serial setting, provides a band matrix preconditioner based on difference
quotients of the backward problem right-hand side function ``fB``. It generates
a banded approximation to the Jacobian with :math:`m_{lB}` sub-diagonals and
:math:`m_{uB}` super-diagonals to be used with one of the Krylov linear solvers.

In order to use the CVBANDPRE module in the solution of the backward problem,
the user need not define any additional functions. Instead, *after* an iterative
``SUNLinearSolver`` object has been attached to CVODES via a call to
:c:func:`CVodeSetLinearSolverB`, the following call to the CVBANDPRE module
initialization function must be made.

.. c:function:: int CVBandPrecInitB(void * cvode_mem, int which, sunindextype nB, sunindextype muB, sunindextype mlB)

   The function :c:func:`CVBandPrecInitB` initializes and allocates  memory for the
   CVBANDPRE preconditioner for the backward problem.  It creates, allocates,
   and stores (internally in the CVODES  solver block) a pointer to the newly
   created CVBANDPRE memory block.

   **Arguments:**
     * ``cvode_mem`` -- pointer to the CVODES memory block.
     * ``which`` -- the identifier of the backward problem.
     * ``nB`` -- backward problem dimension.
     * ``muB`` -- upper half-bandwidth of the backward problem Jacobian approximation.
     * ``mlB`` -- lower half-bandwidth of the backward problem Jacobian approximation.

   **Return value:**
     * ``CVLS_SUCCESS`` -- The call to :c:func:`CVodeBandPrecInitB` was successful.
     * ``CVLS_MEM_FAIL`` -- A memory allocation request has failed.
     * ``CVLS_MEM_NULL`` -- The ``cvode_mem`` argument was ``NULL``.
     * ``CVLS_LMEM_NULL`` -- No linear solver has been attached.
     * ``CVLS_ILL_INPUT`` -- An invalid parameter has been passed.


For more details on CVBANDPRE see :numref:`CVODES.Usage.SIM.precond.cvbandpre`.

Using the band-block-diagonal preconditioner CVBBDPRE
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The adjoint module in CVODES offers an interface to the band-block-diagonal
preconditioner module CVBBDPRE described in section
:numref:`CVODES.Usage.SIM.precond.cvbbdpre`. This generates a preconditioner that
is a block-diagonal matrix with each block being a band matrix and can be used
with one of the Krylov linear solvers and with the MPI-parallel vector module
``NVECTOR_PARALLEL``.

In order to use the CVBBDPRE module in the solution of the backward problem, the
user must define one or two additional functions, described at the end of this
section.

Initialization of CVBBDPRE
""""""""""""""""""""""""""

The CVBBDPRE module is initialized by calling the following function, *after* an
iterative ``SUNLinearSolver`` object has been attached to CVODES via a call to
:c:func:`CVodeSetLinearSolverB`.


.. c:function:: int CVBBDPrecInitB(void * cvode_mem, int which, sunindextype NlocalB, sunindextype mudqB, sunindextype mldqB, sunindextype mukeepB, sunindextype mlkeepB, realtype dqrelyB, CVBBDLocalFnB glocB, CVBBDCommFnB gcommB)

   The function :c:func:`CVBBDPrecInitB` initializes and allocates  memory for the
   CVBBDPRE preconditioner for the backward problem.  It creates, allocates, and
   stores (internally in the CVODES solver  block) a pointer to the newly
   created CVBBDPRE memory block.

   **Arguments:**
     * ``cvode_mem`` -- pointer to the CVODES memory block.
     * ``which`` -- the identifier of the backward problem.
     * ``NlocalB`` -- local vector dimension for the backward problem.
     * ``mudqB`` -- upper half-bandwidth to be used in the difference-quotient Jacobian approximation.
     * ``mldqB`` -- lower half-bandwidth to be used in the difference-quotient Jacobian approximation.
     * ``mukeepB`` -- upper half-bandwidth of the retained banded approximate Jacobian block.
     * ``mlkeepB`` -- lower half-bandwidth of the retained banded approximate Jacobian block.
     * ``dqrelyB`` -- the relative increment in components of ``yB`` used in the difference quotient approximations. The default is :math:`\text{dqrelyB} = \sqrt{\text{unit roundoff}}` , which can be specified by passing ``dqrely = 0.0``.
     * ``glocB`` -- the function which computes the function :math:`g_Bt,y,y_B` approximating the right-hand side of the backward problem.
     * ``gcommB`` -- the optional function which performs all interprocess communication required for the computation of :math:`g_B`.

   **Return value:**
     * ``CVLS_SUCCESS`` -- The call to :c:func:`CVodeBBDPrecInitB` was successful.
     * ``CVLS_MEM_FAIL`` -- A memory allocation request has failed.
     * ``CVLS_MEM_NULL`` -- The ``cvode_mem`` argument was ``NULL``.
     * ``CVLS_LMEM_NULL`` -- No linear solver has been attached.
     * ``CVLS_ILL_INPUT`` -- An invalid parameter has been passed.


.. c:function:: int CVBBDPrecReInitB(void * cvode_mem, int which, sunindextype mudqB, sunindextype mldqB, realtype dqrelyB)

   The function :c:func:`CVBBDPrecReInitB` reinitializes the CVBBDPRE preconditioner
   for the backward problem.

   **Arguments:**
     * ``cvode_mem`` -- pointer to the CVODES memory block returned by :c:func:`CVodeCreate`.
     * ``which`` -- the identifier of the backward problem.
     * ``mudqB`` -- upper half-bandwidth to be used in the difference-quotient Jacobian approximation.
     * ``mldqB`` -- lower half-bandwidth to be used in the difference-quotient Jacobian approximation.
     * ``dqrelyB`` -- the relative increment in components of ``yB`` used in the difference quotient approximations.

   **Return value:**
     * ``CVLS_SUCCESS`` -- The call to :c:func:`CVodeBBDPrecReInitB` was successful.
     * ``CVLS_MEM_FAIL`` -- A memory allocation request has failed.
     * ``CVLS_MEM_NULL`` -- The ``cvode_mem`` argument was ``NULL``.
     * ``CVLS_PMEM_NULL`` -- The :c:func:`CVodeBBDPrecInitB` has not been previously called.
     * ``CVLS_LMEM_NULL`` -- No linear solver has been attached.
     * ``CVLS_ILL_INPUT`` -- An invalid parameter has been passed.


For more details on CVBBDPRE see :numref:`CVODES.Usage.SIM.precond.cvbbdpre`.


User-supplied functions for CVBBDPRE
""""""""""""""""""""""""""""""""""""

To use the CVBBDPRE module, the user must supply one or two functions which the
module calls to construct the preconditioner: a required function ``glocB`` (of
type :c:type:`CVBBDLocalFnB`) which approximates the right-hand side of the
backward problem and which is computed locally, and an optional function
``gcommB`` (of type :c:type:`CVBBDCommFnB`) which performs all interprocess
communication necessary to evaluate this approximate right-hand side. The
prototypes for these two functions are described below.


.. c:type:: int (*CVBBDLocalFnB)(sunindextype NlocalB, realtype t, N_Vector y, N_Vector yB, N_Vector gB, void *user_dataB)

   This ``glocB`` function loads the vector ``gB``, an approximation to the
   right-hand side :math:`f_B` of the backward problem, as a function of ``t``,
   ``y``,  and ``yB``.

   **Arguments:**
     * ``NlocalB`` -- is the local vector length for the backward problem.
     * ``t`` -- is the value of the independent variable.
     * ``y`` -- is the current value of the forward solution vector.
     * ``yB`` -- is the current value of the backward dependent variable vector.
     * ``gB`` -- is the output vector, :math:`g_B(t, y, y_B)`.
     * ``user_dataB`` -- is a pointer to the same user data passed to :c:func:`CVodeSetUserDataB`.

   **Return value:**
      An :c:type:`CVBBDLocalFnB` should return 0 if successful, a positive value
      if a recoverable error occurred (in which case CVODES will attempt to
      correct), or a negative value if it failed unrecoverably (in which case
      the integration is halted and :c:func:`CVodeB` returns
      ``CV_LSETUP_FAIL``).

   **Notes:**
      This routine must assume that all interprocess communication of data
      needed to  calculate ``gB`` has already been done, and this data is
      accessible within  ``user_dataB``.

   .. warning::
      Before calling the user's :c:type:`CVBBDLocalFnB`, CVODES needs to
      evaluate  (through interpolation) the values of the states from the
      forward integration.  If an error occurs in the interpolation, CVODES
      triggers an unrecoverable  failure in the preconditioner setup function
      which will halt the integration  (:c:func:`CVodeB` returns
      ``CV_LSETUP_FAIL``).


.. c:type:: int (*CVBBDCommFnB)(sunindextype NlocalB, realtype t, N_Vector y, N_Vector yB, void *user_dataB)

   This ``gcommB`` function must perform all interprocess communications
   necessary  for the execution of the ``glocB`` function above, using the input
   vectors ``y`` and ``yB``.

   **Arguments:**
     * ``NlocalB`` -- is the local vector length.
     * ``t`` -- is the value of the independent variable.
     * ``y`` -- is the current value of the forward solution vector.
     * ``yB`` -- is the current value of the backward dependent variable vector.
     * ``user_dataB`` -- is a pointer to the same user data passed to :c:func:`CVodeSetUserDataB`.

   **Return value:**
      An :c:type:`CVBBDCommFnB` should return 0 if successful, a positive value
      if a recoverable error occurred (in which case CVODES will attempt to
      correct), or a negative value if it failed unrecoverably (in which case
      the integration is halted and :c:func:`CVodeB` returns
      ``CV_LSETUP_FAIL``).

   **Notes:**
      The ``gcommB`` function is expected to save communicated data in space
      defined within the  structure ``user_dataB``.  Each call to the ``gcommB``
      function is preceded by a call to the function that  evaluates the
      right-hand side of the backward problem with the same ``t``, ``y``,  and
      ``yB``, arguments. If there is no additional communication needed, then
      pass ``gcommB = NULL`` to :c:func:`CVBBDPrecInitB`.
