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

.. _ARKODE.Usage.ARKStep.PreconditionerModules:

Preconditioner modules
============================

The efficiency of Krylov iterative methods for the solution of linear
systems can be greatly enhanced through preconditioning.  For problems
in which the user cannot define a more effective, problem-specific
preconditioner, ARKODE provides two internal preconditioner modules
that may be used by ARKStep: a banded preconditioner for serial and
threaded problems (ARKBANDPRE) and a band-block-diagonal
preconditioner for parallel problems (ARKBBDPRE).


.. _ARKODE.Usage.ARKStep.BandPre:

A serial banded preconditioner module
-------------------------------------------

This preconditioner provides a band matrix preconditioner for use with
iterative SUNLINSOL modules in a serial or threaded setting. It requires
that the problem be set up using either the
NVECTOR_SERIAL, NVECTOR_OPENMP or NVECTOR_PTHREADS module, due to data
access patterns.  It also currently requires that the problem involve
an identity mass matrix, i.e., :math:`M = I`.

This module uses difference quotients of the ODE right-hand
side function :math:`f^I` to generate a band matrix of bandwidth
``ml + mu + 1``, where the number of super-diagonals (``mu``, the
upper half-bandwidth) and sub-diagonals (``ml``, the lower
half-bandwidth) are specified by the user.  This band matrix is used
to to form a preconditioner the Krylov linear solver.  Although this
matrix is intended to approximate the Jacobian
:math:`J = \dfrac{\partial f^I}{\partial y}`, it may be a very crude
approximation, since the true Jacobian may not be banded, or its true
bandwidth may be larger than ``ml + mu + 1``.  However, as long as the
banded approximation generated for the preconditioner is sufficiently
accurate, it may speed convergence of the Krylov iteration.



ARKBANDPRE usage
"""""""""""""""""""""

In order to use the ARKBANDPRE module, the user need not define
any additional functions.  In addition to the header files required
for the integration of the ODE problem (see
:numref:`ARKODE.Usage.Headers`), to use the ARKBANDPRE module, the user's
program must include the header file ``arkode_bandpre.h`` which
declares the needed function prototypes.  The following is a summary
of the usage of this module.  Steps that are unchanged from the
skeleton program presented in :numref:`ARKODE.Usage.ARKStep.Skeleton` are
*italicized*.

1. *Initialize multi-threaded environment (if appropriate)*

2. *Set problem dimensions*

3. *Set vector of initial values*

4. *Create ARKStep object*

5. *Specify integration tolerances*

6. Create iterative linear solver object

   When creating the iterative linear solver object, specify the type
   of preconditioning (``SUN_PREC_LEFT`` or ``SUN_PREC_RIGHT``) to use.

7. *Set linear solver optional inputs*

8. *Attach linear solver module*

9. Initialize the ARKBANDPRE preconditioner module

    Specify the upper and lower half-bandwidths (``mu`` and ``ml``,
    respectively) and call

    ``ier = ARKBandPrecInit(arkode_mem, N, mu, ml);``

    to allocate memory and initialize the internal preconditioner
    data.

10. *Set optional inputs*

    Note that the user should not call
    :c:func:`ARKStepSetPreconditioner()` as it will overwrite the
    preconditioner setup and solve functions.

11. *Create nonlinear solver object*

12. *Attach nonlinear solver module*

13. *Set nonlinear solver optional inputs*

14. *Specify rootfinding problem*

15. *Advance solution in time*

16. Get optional outputs

    Additional optional outputs associated with ARKBANDPRE are
    available by way of the two routines described below,
    :c:func:`ARKBandPrecGetWorkSpace()` and
    :c:func:`ARKBandPrecGetNumRhsEvals()`.

17. *Deallocate memory for solution vector*

18. *Free solver memory*

19. *Free linear solver memory*





ARKBANDPRE user-callable functions
"""""""""""""""""""""""""""""""""""""

The ARKBANDPRE preconditioner module is initialized and attached
by calling the following function:



.. c:function:: int ARKBandPrecInit(void* arkode_mem, sunindextype N, sunindextype mu, sunindextype ml)

   Initializes the ARKBANDPRE preconditioner and
   allocates required (internal) memory for it.

   **Arguments:**
      * *arkode_mem* -- pointer to the ARKStep memory block.
      * *N* -- problem dimension (size of ODE system).
      * *mu* -- upper half-bandwidth of the Jacobian approximation.
      * *ml* -- lower half-bandwidth of the Jacobian approximation.

   **Return value:**
      * *ARKLS_SUCCESS* if no errors occurred
      * *ARKLS_MEM_NULL* if the ARKStep memory is ``NULL``
      * *ARKLS_LMEM_NULL* if the linear solver memory is ``NULL``
      * *ARKLS_ILL_INPUT* if an input has an illegal value
      * *ARKLS_MEM_FAIL* if a memory allocation request failed

   **Notes:**
      The banded approximate Jacobian will have nonzero elements
      only in locations :math:`(i,j)` with *ml* :math:`\le j-i \le` *mu*.



The following two optional output functions are available for use with
the ARKBANDPRE module:



.. c:function:: int ARKBandPrecGetWorkSpace(void* arkode_mem, long int* lenrwLS, long int* leniwLS)

   Returns the sizes of the ARKBANDPRE real and integer
   workspaces.

   **Arguments:**
      * *arkode_mem* -- pointer to the ARKStep memory block.
      * *lenrwLS* -- the number of ``realtype`` values in the
        ARKBANDPRE workspace.
      * *leniwLS* -- the number of integer values in the  ARKBANDPRE workspace.

   **Return value:**
      * *ARKLS_SUCCESS* if no errors occurred
      * *ARKLS_MEM_NULL* if the ARKStep memory is ``NULL``
      * *ARKLS_LMEM_NULL* if the linear solver memory is ``NULL``
      * *ARKLS_PMEM_NULL* if the preconditioner memory is ``NULL``

   **Notes:**
      The workspace requirements reported by this routine
      correspond only to memory allocated within the ARKBANDPRE module
      (the banded matrix approximation, banded ``SUNLinearSolver``
      object, and temporary vectors).

      The workspaces referred to here exist in addition to those given by
      the corresponding function :c:func:`ARKStepGetLSWorkspace()`.



.. c:function:: int ARKBandPrecGetNumRhsEvals(void* arkode_mem, long int* nfevalsBP)

   Returns the number of calls made to the user-supplied
   right-hand side function :math:`f^I` for constructing the
   finite-difference banded Jacobian approximation used within the
   preconditioner setup function.

   **Arguments:**
      * *arkode_mem* -- pointer to the ARKStep memory block.
      * *nfevalsBP* -- number of calls to :math:`f^I`.

   **Return value:**
      * *ARKLS_SUCCESS* if no errors occurred
      * *ARKLS_MEM_NULL* if the ARKStep memory is ``NULL``
      * *ARKLS_LMEM_NULL* if the linear solver memory is ``NULL``
      * *ARKLS_PMEM_NULL* if the preconditioner memory is ``NULL``

   **Notes:**
      The counter *nfevalsBP* is distinct from the counter
      *nfevalsLS* returned by the corresponding function
      :c:func:`ARKStepGetNumLSRhsEvals()` and also from *nfi_evals* returned by
      :c:func:`ARKStepGetNumRhsEvals()`.  The total number of right-hand
      side function evaluations is the sum of all three of these
      counters, plus the *nfe_evals* counter for :math:`f^E` calls
      returned by :c:func:`ARKStepGetNumRhsEvals()`.





.. _ARKODE.Usage.ARKStep.BBDPre:

A parallel band-block-diagonal preconditioner module
---------------------------------------------------------

A principal reason for using a parallel ODE solver (such as ARKODE)
lies in the solution of partial differential equations
(PDEs). Moreover, Krylov iterative methods are used on many such
problems due to the nature of the underlying linear system of
equations that needs to solved at each time step.  For many PDEs, the
linear algebraic system is large, sparse and structured.  However, if
a Krylov iterative method is to be effective in this setting, then a
nontrivial preconditioner is required.  Otherwise, the rate of
convergence of the Krylov iterative method is usually slow, and
degrades as the PDE mesh is refined.  Typically, an effective
preconditioner must be problem-specific.

However, we have developed one type of preconditioner that treats a
rather broad class of PDE-based problems.  It has been successfully
used with CVODE for several realistic, large-scale problems :cite:p:`HiTa:98`.
It is included in a software module within the ARKODE package, and is
accessible within the ARKStep time stepping module.  This
preconditioning module works with the parallel vector module
NVECTOR_PARALLEL and is usable with any of the Krylov iterative linear
solvers through the ARKLS interface. It generates a preconditioner
that is a block-diagonal matrix with each block being a band
matrix. The blocks need not have the same number of super- and
sub-diagonals and these numbers may vary from block to block. This
Band-Block-Diagonal Preconditioner module is called ARKBBDPRE.

One way to envision these preconditioners is to think of the
computational PDE domain as being subdivided into :math:`Q`
non-overlapping subdomains, where each subdomain is assigned to one of
the :math:`Q` MPI tasks used to solve the ODE system.  The basic idea
is to isolate the preconditioning so that it is local to each process,
and also to use a (possibly cheaper) approximate right-hand side
function for construction of this preconditioning matrix.  This
requires the definition of a new function :math:`g(t,y) \approx
f^I(t,y)` that will be used to construct the BBD preconditioner
matrix.  At present, we assume that the ODE be written in explicit
form as

.. math::
   \dot{y} = f^E(t,y) + f^I(t,y),

where :math:`f^I` corresponds to the ODE components to be treated
implicitly, i.e. this preconditioning module does not support problems
with non-identity mass matrices.  The user may set :math:`g = f^I`, if
no less expensive approximation is desired.

Corresponding to the domain decomposition, there is a decomposition of
the solution vector :math:`y` into :math:`Q` disjoint blocks
:math:`y_q`, and a decomposition of :math:`g` into blocks
:math:`g_q`. The block :math:`g_q` depends both on :math:`y_p` and on
components of blocks :math:`y_{q'}` associated with neighboring
subdomains (so-called ghost-cell data).  If we let :math:`\bar{y}_q`
denote :math:`y_q` augmented with those other components on which
:math:`g_q` depends, then we have

.. math::
   g(t,y) = \left[ g_1(t,\bar{y}_1), g_2(t,\bar{y}_2), \ldots , g_Q(t,\bar{y}_Q) \right]^T,

and each of the blocks :math:`g_q(t,\bar{y}_q)` is decoupled from one another.

The preconditioner associated with this decomposition has the form

.. math::
   P = \begin{bmatrix} P_1 & & & \\ & P_2 & & \\ & & \ddots &\\ & & & P_Q \end{bmatrix}

where

.. math::
   P_q \approx I - \gamma J_q

and where :math:`J_q` is a difference quotient approximation to
:math:`\dfrac{\partial g_q}{\partial \bar{y}_q}`.  This matrix is taken
to be banded, with upper and lower half-bandwidths *mudq* and
*mldq* defined as the number of non-zero diagonals above and below
the main diagonal, respectively.  The difference quotient
approximation is computed using *mudq* + *mldq* + 2 evaluations of
:math:`g_m`, but only a matrix of bandwidth *mukeep* + *mlkeep* + 1 is
retained. Neither pair of parameters need be the true half-bandwidths
of the Jacobian of the local block of :math:`g`, if smaller values
provide a more efficient preconditioner. The solution of the complete
linear system

.. math::
   Px = b

reduces to solving each of the distinct equations

.. math::
   P_q x_q = b_q, \quad q=1,\ldots,Q,

and this is done by banded LU factorization of :math:`P_q` followed by
a banded backsolve.

Similar block-diagonal preconditioners could be considered with
different treatments of the blocks :math:`P_q`.  For example,
incomplete LU factorization or an iterative method could be used
instead of banded LU factorization.



ARKBBDPRE user-supplied functions
""""""""""""""""""""""""""""""""""

The ARKBBDPRE module calls two user-provided functions to construct
:math:`P`: a required function *gloc* (of type :c:func:`ARKLocalFn()`)
which approximates the right-hand side function :math:`g(t,y) \approx
f^I(t,y)` and which is computed locally, and an optional function
*cfn* (of type :c:func:`ARKCommFn()`) which performs all inter-process
communication necessary to evaluate the approximate right-hand side
:math:`g`. These are in addition to the user-supplied right-hand side
function :math:`f^I`. Both functions take as input the same pointer
*user_data* that is passed by the user to
:c:func:`ARKStepSetUserData()` and that was passed to the user's
function :math:`f^I`. The user is responsible for providing space
(presumably within *user_data*) for components of :math:`y` that are
communicated between processes by *cfn*, and that are then used by
*gloc*, which should not do any communication.



.. c:type:: int (*ARKLocalFn)(sunindextype Nlocal, realtype t, N_Vector y, N_Vector glocal, void* user_data)

   This *gloc* function computes :math:`g(t,y)`.  It
   fills the vector *glocal* as a function of *t* and *y*.

   **Arguments:**
      * *Nlocal* -- the local vector length.
      * *t* -- the value of the independent variable.
      * *y* -- the value of the dependent variable vector on this process.
      * *glocal* -- the output vector of :math:`g(t,y)` on this process.
      * *user_data* -- a pointer to user data, the same as the
        *user_data* parameter passed to :c:func:`ARKStepSetUserData()`.

   **Return value:**
      An *ARKLocalFn* should return 0 if successful, a positive value if
      a recoverable error occurred (in which case ARKStep will attempt to
      correct), or a negative value if it failed unrecoverably (in which
      case the integration is halted and :c:func:`ARKStepEvolve()` will return
      *ARK_LSETUP_FAIL*).

   **Notes:**
      This function should assume that all inter-process
      communication of data needed to calculate *glocal* has already been
      done, and that this data is accessible within user data.

      The case where :math:`g` is mathematically identical to :math:`f^I`
      is allowed.



.. c:type:: int (*ARKCommFn)(sunindextype Nlocal, realtype t, N_Vector y, void* user_data)

   This *cfn* function performs all inter-process
   communication necessary for the execution of the *gloc* function
   above, using the input vector *y*.

   **Arguments:**
      * *Nlocal* -- the local vector length.
      * *t* -- the value of the independent variable.
      * *y* -- the value of the dependent variable vector on this process.
      * *user_data* -- a pointer to user data, the same as the
        *user_data* parameter passed to :c:func:`ARKStepSetUserData()`.

   **Return value:**
      An *ARKCommFn* should return 0 if successful, a positive value if a
      recoverable error occurred (in which case ARKStep will attempt to
      correct), or a negative value if it failed unrecoverably (in which
      case the integration is halted and :c:func:`ARKStepEvolve()` will return
      *ARK_LSETUP_FAIL*).

   **Notes:**
      The *cfn* function is expected to save communicated data in
      space defined within the data structure *user_data*.

      Each call to the *cfn* function is preceded by a call to the
      right-hand side function :math:`f^I` with the same :math:`(t,y)`
      arguments. Thus, *cfn* can omit any communication done by
      :math:`f^I` if relevant to the evaluation of *glocal*. If all
      necessary communication was done in :math:`f^I`, then *cfn* =
      ``NULL`` can be passed in the call to :c:func:`ARKBBDPrecInit()`
      (see below).




ARKBBDPRE usage
"""""""""""""""""""""

In addition to the header files required for the integration of the
ODE problem (see :numref:`ARKODE.Usage.Headers`), to use the
ARKBBDPRE module, the user's program must include the header file
``arkode_bbdpre.h`` which declares the needed function prototypes.

The following is a summary of the proper usage of this module. Steps
that are unchanged from the skeleton program presented in
:numref:`ARKODE.Usage.ARKStep.Skeleton` are *italicized*.

1. *Initialize MPI*

2. *Set problem dimensions*

3. *Set vector of initial values*

4. *Create ARKStep object*

5. *Specify integration tolerances*

6. Create iterative linear solver object

   When creating the iterative linear solver object, specify the type
   of preconditioning (``SUN_PREC_LEFT`` or ``SUN_PREC_RIGHT``) to use.

7. *Set linear solver optional inputs*

8. *Attach linear solver module*

9. Initialize the ARKBBDPRE preconditioner module

   Specify the upper and lower half-bandwidths for computation
   ``mudq`` and ``mldq``, the upper and lower half-bandwidths for
   storage ``mukeep`` and ``mlkeep``, and call

   ``ier = ARKBBDPrecInit(arkode_mem, Nlocal, mudq, mldq, mukeep, mlkeep, dqrely, gloc, cfn);``

   to allocate memory and initialize the internal preconditioner
   data. The last two arguments of :c:func:`ARKBBDPrecInit()` are the
   two user-supplied functions of type :c:func:`ARKLocalFn()` and
   :c:func:`ARKCommFn()` described above, respectively.

10. *Set optional inputs*

    Note that the user should not call
    :c:func:`ARKStepSetPreconditioner()` as it will overwrite the
    preconditioner setup and solve functions.

11. *Create nonlinear solver object*

12. *Attach nonlinear solver module*

13. *Set nonlinear solver optional inputs*

14. *Specify rootfinding problem*

15. *Advance solution in time*

16. *Get optional outputs*

    Additional optional outputs associated with ARKBBDPRE are
    available through the routines
    :c:func:`ARKBBDPrecGetWorkSpace()` and
    :c:func:`ARKBBDPrecGetNumGfnEvals()`.

17. *Deallocate memory for solution vector*

18. *Free solver memory*

19. *Free linear solver memory*

20. *Finalize MPI*





ARKBBDPRE user-callable functions
""""""""""""""""""""""""""""""""""""

The ARKBBDPRE preconditioner module is initialized (or re-initialized)
and attached to the integrator by calling the following functions:

.. c:function:: int ARKBBDPrecInit(void* arkode_mem, sunindextype Nlocal, sunindextype mudq, sunindextype mldq, sunindextype mukeep, sunindextype mlkeep, realtype dqrely, ARKLocalFn gloc, ARKCommFn cfn)

   Initializes and allocates (internal) memory for the
   ARKBBDPRE preconditioner.

   **Arguments:**
      * *arkode_mem* -- pointer to the ARKStep memory block.
      * *Nlocal* -- local vector length.
      * *mudq* -- upper half-bandwidth to be used in the difference
        quotient Jacobian approximation.
      * *mldq* -- lower half-bandwidth to be used in the difference
        quotient Jacobian approximation.
      * *mukeep* -- upper half-bandwidth of the retained banded
        approximate Jacobian block.
      * *mlkeep* -- lower half-bandwidth of the retained banded
        approximate Jacobian block.
      * *dqrely* -- the relative increment in components of *y* used in
        the difference quotient approximations.  The default is *dqrely*
        = :math:`\sqrt{\text{unit roundoff}}`, which can be specified by
        passing *dqrely* = 0.0.
      * *gloc* -- the name of the C function (of type :c:func:`ARKLocalFn()`)
        which computes the approximation :math:`g(t,y) \approx f^I(t,y)`.
      * *cfn* -- the name of the C function (of type :c:func:`ARKCommFn()`) which
        performs all inter-process communication required for the
        computation of :math:`g(t,y)`.

   **Return value:**
      * *ARKLS_SUCCESS* if no errors occurred
      * *ARKLS_MEM_NULL* if the ARKStep memory is ``NULL``
      * *ARKLS_LMEM_NULL* if the linear solver memory is ``NULL``
      * *ARKLS_ILL_INPUT* if an input has an illegal value
      * *ARKLS_MEM_FAIL* if a memory allocation request failed

   **Notes:**
      If one of the half-bandwidths *mudq* or *mldq* to be used
      in the difference quotient calculation of the approximate Jacobian is
      negative or exceeds the value *Nlocal*-1, it is replaced by 0 or
      *Nlocal*-1 accordingly.

      The half-bandwidths *mudq* and *mldq* need not be the true
      half-bandwidths of the Jacobian of the local block of :math:`g`
      when smaller values may provide a greater efficiency.

      Also, the half-bandwidths *mukeep* and *mlkeep* of the retained
      banded approximate Jacobian block may be even smaller than
      *mudq* and *mldq*, to reduce storage and computational costs
      further.

      For all four half-bandwidths, the values need not be the same on
      every processor.



The ARKBBDPRE module also provides a re-initialization function to
allow solving a sequence of problems of the same size, with the same
linear solver choice, provided there is no change in *Nlocal*,
*mukeep*, or *mlkeep*. After solving one problem, and after
calling :c:func:`ARKStepReInit()` to re-initialize ARKStep for a
subsequent problem, a call to :c:func:`ARKBBDPrecReInit()` can be made
to change any of the following: the half-bandwidths *mudq* and
*mldq* used in the difference-quotient Jacobian approximations, the
relative increment *dqrely*, or one of the user-supplied functions
*gloc* and *cfn*. If there is a change in any of the linear solver
inputs, an additional call to the "Set" routines provided by the
SUNLINSOL module, and/or one or more of the corresponding
``ARKStepSet***`` functions, must also be made (in the proper order).


.. c:function:: int ARKBBDPrecReInit(void* arkode_mem, sunindextype mudq, sunindextype mldq, realtype dqrely)

   Re-initializes the ARKBBDPRE preconditioner module.

   **Arguments:**
      * *arkode_mem* -- pointer to the ARKStep memory block.
      * *mudq* -- upper half-bandwidth to be used in the difference
        quotient Jacobian approximation.
      * *mldq* -- lower half-bandwidth to be used in the difference
        quotient Jacobian approximation.
      * *dqrely* -- the relative increment in components of *y* used in
        the difference quotient approximations.  The default is *dqrely*
        = :math:`\sqrt{\text{unit roundoff}}`, which can be specified by
        passing *dqrely* = 0.0.

   **Return value:**
      * *ARKLS_SUCCESS* if no errors occurred
      * *ARKLS_MEM_NULL* if the ARKStep memory is ``NULL``
      * *ARKLS_LMEM_NULL* if the linear solver memory is ``NULL``
      * *ARKLS_PMEM_NULL* if the preconditioner memory is ``NULL``

   **Notes:**
      If one of the half-bandwidths *mudq* or *mldq* is
      negative or exceeds the value *Nlocal*-1, it is replaced by 0 or
      *Nlocal*-1 accordingly.


The following two optional output functions are available for use with
the ARKBBDPRE module:


.. c:function:: int ARKBBDPrecGetWorkSpace(void* arkode_mem, long int* lenrwBBDP, long int* leniwBBDP)

   Returns the processor-local ARKBBDPRE real and
   integer workspace sizes.

   **Arguments:**
      * *arkode_mem* -- pointer to the ARKStep memory block.
      * *lenrwBBDP* -- the number of ``realtype`` values in the
        ARKBBDPRE workspace.
      * *leniwBBDP* -- the number of integer values in the  ARKBBDPRE workspace.

   **Return value:**
      * *ARKLS_SUCCESS* if no errors occurred
      * *ARKLS_MEM_NULL* if the ARKStep memory is ``NULL``
      * *ARKLS_LMEM_NULL* if the linear solver memory is ``NULL``
      * *ARKLS_PMEM_NULL* if the preconditioner memory is ``NULL``

   **Notes:**
      The workspace requirements reported by this routine
      correspond only to memory allocated within the ARKBBDPRE module
      (the banded matrix approximation, banded ``SUNLinearSolver``
      object, temporary vectors). These values are local to each process.

      The workspaces referred to here exist in addition to those given by
      the corresponding function :c:func:`ARKStepGetLSWorkSpace()`.



.. c:function:: int ARKBBDPrecGetNumGfnEvals(void* arkode_mem, long int* ngevalsBBDP)

   Returns the number of calls made to the user-supplied
   *gloc* function (of type :c:func:`ARKLocalFn()`) due to the finite
   difference approximation of the Jacobian blocks used within the
   preconditioner setup function.

   **Arguments:**
      * *arkode_mem* -- pointer to the ARKStep memory block.
      * *ngevalsBBDP* -- the number of calls made to the user-supplied
        *gloc* function.

   **Return value:**
      * *ARKLS_SUCCESS* if no errors occurred
      * *ARKLS_MEM_NULL* if the ARKStep memory is ``NULL``
      * *ARKLS_LMEM_NULL* if the linear solver memory is ``NULL``
      * *ARKLS_PMEM_NULL* if the preconditioner memory is ``NULL``


In addition to the *ngevalsBBDP* *gloc* evaluations, the costs
associated with ARKBBDPRE also include *nlinsetups* LU
factorizations, *nlinsetups* calls to *cfn*, *npsolves* banded
backsolve calls, and *nfevalsLS* right-hand side function
evaluations, where *nlinsetups* is an optional ARKStep output and
*npsolves* and *nfevalsLS* are linear solver optional outputs (see
the table :numref:`ARKODE.Usage.ARKStep.ARKLsOutputs`).
