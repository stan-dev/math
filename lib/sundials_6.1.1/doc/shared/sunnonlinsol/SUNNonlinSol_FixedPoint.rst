..
   Programmer(s): Daniel R. Reynolds @ SMU
   ----------------------------------------------------------------
   SUNDIALS Copyright Start
   Copyright (c) 2002-2022, Lawrence Livermore National Security
   and Southern Methodist University.
   All rights reserved.

   See the top-level LICENSE and NOTICE files for details.

   SPDX-License-Identifier: BSD-3-Clause
   SUNDIALS Copyright End
   ----------------------------------------------------------------

.. _SUNNonlinSol.FixedPoint:

==================================================
The SUNNonlinSol_FixedPoint implementation
==================================================

This section describes the SUNNonlinSol implementation of a fixed point
(functional) iteration with optional Anderson acceleration. To access the
SUNNonlinSol_FixedPoint module, include the header file
``sunnonlinsol/sunnonlinsol_fixedpoint.h``. We note that the
SUNNonlinSol_FixedPoint module is accessible from SUNDIALS integrators
*without* separately linking to the
``libsundials_sunnonlinsolfixedpoint`` module library.


.. _SUNNonlinSol.FixedPoint.Math:

SUNNonlinSol_FixedPoint description
-----------------------------------------------

To find the solution to

.. math::
   G(y) = y \,
   :label: e:fixed_point_sys

given an initial guess :math:`y^{(0)}`, the fixed point iteration
computes a series of approximate solutions

.. math::
   y^{(n+1)} = G(y^{(n)})
   :label: e:fixed_point_iteration

where :math:`n` is the iteration index. The convergence of this
iteration may be accelerated using Anderson's method
:cite:p:`Anderson65,Walker-Ni09,Fang-Saad09,LWWY11`.  With Anderson
acceleration using subspace size :math:`m`, the series of approximate
solutions can be formulated as the linear combination

.. math::
   y^{(n+1)} = \beta \sum_{i=0}^{m_n} \alpha_i^{(n)} G(y^{(n-m_n+i)}) + (1 - \beta) \sum_{i=0}^{m_n} \alpha_i^{(n)} y_{n-m_n+i}
   :label: e:accelerated_fixed_point_iteration

where :math:`m_n = \min{\{m,n\}}` and the factors

.. math::
   \alpha^{(n)} =(\alpha_0^{(n)}, \ldots, \alpha_{m_n}^{(n)})

solve the minimization problem :math:`\min\limits_\alpha  \| F_n \alpha^T
\|_2` under the constraint that :math:`\sum\limits_{i=0}^{m_n} \alpha_i = 1` where

.. math::
   F_{n} = (f_{n-m_n}, \ldots, f_{n})

with :math:`f_i = G(y^{(i)}) - y^{(i)}`. Due to this constraint, in
the limit of :math:`m=0` the accelerated fixed point iteration formula
:eq:`e:accelerated_fixed_point_iteration` simplifies to the standard
fixed point iteration :eq:`e:fixed_point_iteration`.

Following the recommendations made in :cite:p:`Walker-Ni09`, the
SUNNonlinSol_FixedPoint implementation computes the series of
approximate solutions as

.. math::
   y^{(n+1)} = G(y^{(n)})-\sum_{i=0}^{m_n-1} \gamma_i^{(n)} \Delta g_{n-m_n+i} - (1 - \beta) (f(y^{(n)}) - \sum_{i=0}^{m_n-1} \gamma_i^{(n)} \Delta f_{n-m_n+i})
   :label: e:accelerated_fixed_point_iteration_impl

with :math:`\Delta g_i = G(y^{(i+1)}) - G(y^{(i)})` and where the
factors

.. math::
   \gamma^{(n)} =(\gamma_0^{(n)}, \ldots, \gamma_{m_n-1}^{(n)})

solve the unconstrained minimization problem
:math:`\min\limits_\gamma \| f_n - \Delta F_n \gamma^T \|_2` where

.. math::
   \Delta F_{n} = (\Delta f_{n-m_n}, \ldots, \Delta f_{n-1}),

with :math:`\Delta f_i = f_{i+1} - f_i`. The least-squares problem is
solved by applying a QR factorization to :math:`\Delta F_n = Q_n R_n`
and solving  :math:`R_n \gamma = Q_n^T f_n`.

The acceleration subspace size :math:`m` is required when constructing
the SUNNonlinSol_FixedPoint object.  The default maximum number of
iterations and the stopping criteria for the fixed point iteration are
supplied by the SUNDIALS integrator when SUNNonlinSol_FixedPoint
is attached to it.  Both the maximum number of iterations and the
convergence test function may be modified by the user by calling
:c:func:`SUNNonlinSolSetMaxIters` and
:c:func:`SUNNonlinSolSetConvTestFn` after attaching the
SUNNonlinSol_FixedPoint object to the integrator.


.. _SUNNonlinSol.FixedPoint.Functions:

SUNNonlinSol_FixedPoint functions
--------------------------------------------

The SUNNonlinSol_FixedPoint module provides the following constructor
for creating the ``SUNNonlinearSolver`` object.



.. c:function:: SUNNonlinearSolver SUNNonlinSol_FixedPoint(N_Vector y, int m, SUNContext sunctx)

   This creates a ``SUNNonlinearSolver`` object for use with SUNDIALS
   integrators to solve nonlinear systems of the form :math:`G(y) = y`.

   **Arguments:**
      * *y* -- a template for cloning vectors needed within the solver.
      * *m* -- the number of acceleration vectors to use.
      * *sunctx* -- the :c:type:`SUNContext` object (see :numref:`SUNDIALS.SUNContext`)

   **Return value:**
      A SUNNonlinSol object if the constructor exits successfully,
      otherwise it will be ``NULL``.


Since the accelerated fixed point iteration
:eq:`e:fixed_point_iteration` does not require the setup or solution
of any linear systems, the SUNNonlinSol_FixedPoint module implements
all of the functions defined in
:numref:`SUNNonlinSol.API.CoreFn`--:numref:`SUNNonlinSol.API.GetFn`
except for the :c:func:`SUNNonlinSolSetup`,
:c:func:`SUNNonlinSolSetLSetupFn`, and :c:func:`SUNNonlinSolSetLSolveFn`
functions, that are set to ``NULL``. The SUNNonlinSol_FixedPoint
functions have the same names as those defined by the generic
SUNNonlinSol API with ``_FixedPoint`` appended to the function name.
Unless using the SUNNonlinSol_FixedPoint module as a standalone
nonlinear solver the generic functions defined in
:numref:`SUNNonlinSol.API.CoreFn`--:numref:`SUNNonlinSol.API.GetFn`
should be called in favor of the SUNNonlinSol_FixedPoint-specific
implementations.

The SUNNonlinSol_FixedPoint module also defines the following
user-callable functions.



.. c:function:: int SUNNonlinSolGetSysFn_FixedPoint(SUNNonlinearSolver NLS, SUNNonlinSolSysFn *SysFn)

   This returns the fixed-point function that defines the nonlinear system.

   **Arguments:**
      * *NLS* -- a SUNNonlinSol object.
      * *SysFn* -- the function defining the nonlinear system.

   **Return value:**
      The return value is zero for a successful call, and a
      negative value for a failure.

   **Notes:**
      This function is intended for users that wish to evaluate the
      fixed-point function in a custom convergence test function for
      the SUNNonlinSol_FixedPoint module. We note that
      SUNNonlinSol_FixedPoint will not leverage the results from any user
      calls to *SysFn*.


.. c:function:: int SUNNonlinSolSetDamping_FixedPoint(SUNNonlinearSolver NLS, realtype beta)

   This sets the damping parameter :math:`\beta` to use with Anderson
   acceleration. By default damping is disabled i.e., :math:`\beta = 1.0`.

   **Arguments:**
     * *NLS* -- a SUNNonlinSol object.
     * *beta* -- the damping parameter :math:`0 < \beta \leq 1`.

   **Return value:**
      * ``SUN_NLS_SUCCESS`` if successful.
      * ``SUN_NLS_MEM_NULL`` if ``NLS`` was ``NULL``.
      * ``SUN_NLS_ILL_INPUT`` if ``beta`` was negative.

   **Notes:**
      A ``beta`` value should satisfy :math:`0 < \beta < 1` if
      damping is to be used. A value of one or more will disable damping.


.. c:function:: int SUNNonlinSolSetInfoFile_FixedPoint(SUNNonlinearSolver NLS, FILE* info_file)

   Thissets the output file where all informative (non-error)
   messages should be directed.

   **Arguments:**
      * *NLS* -- a SUNNonlinSol object.
      * *info_file* -- pointer to output file (``stdout`` by default);
         a ``NULL`` input will disable output.

   **Return value:**
      * ``SUN_NLS_SUCCESS`` if successful.
      * ``SUN_NLS_MEM_NULL`` if ``NLS`` was ``NULL``.
      * ``SUN_NLS_ILL_INPUT`` if SUNDIALS was not built with monitoring enabled.

   **Notes:**
      This function is intended for users that wish to monitor the nonlinear
      solver progress. By default, the file pointer is set to ``stdout``.

   .. warning::

      SUNDIALS must be built with the CMake option
      ``SUNDIALS_BUILD_WITH_MONITORING`` to utilize this function.
      See :numref:`Installation.CMake.Options` for more information.


.. c:function:: int SUNNonlinSolSetPrintLevel_FixedPoint(SUNNonlinearSolver NLS, int print_level)

   This specifies the level of verbosity of the output.

   **Arguments:**
      * *NLS* -- a SUNNonlinSol object.
      * *print_level* -- flag indicating level of verbosity;
        must be one of:

         * 0, no information is printed (default).
         * 1, for each nonlinear iteration the residual norm is printed.

   **Return value:**
      * ``SUN_NLS_SUCCESS`` if successful.
      * ``SUN_NLS_MEM_NULL`` if ``NLS`` was ``NULL``.
      * ``SUN_NLS_ILL_INPUT`` if SUNDIALS was not built with monitoring enabled,
        or the print level value was invalid.

   **Notes:**
      This function is intended for users that wish to monitor the nonlinear
      solver progress. By default, the print level is 0.

   .. warning::

      SUNDIALS must be built with the CMake option
      ``SUNDIALS_BUILD_WITH_MONITORING`` to utilize this function.
      See :numref:`Installation.CMake.Options` for more information.


.. _SUNNonlinSol.FixedPoint.Content:

SUNNonlinSol_FixedPoint content
----------------------------------------

The *content* field of the SUNNonlinSol_FixedPoint module is the
following structure.

.. code-block:: c

   struct _SUNNonlinearSolverContent_FixedPoint {

     SUNNonlinSolSysFn      Sys;
     SUNNonlinSolConvTestFn CTest;

     int          m;
     int         *imap;
     realtype    *R;
     booleantype  damping
     realtype     beta
     realtype    *gamma;
     realtype    *cvals;
     N_Vector    *df;
     N_Vector    *dg;
     N_Vector    *q;
     N_Vector    *Xvecs;
     N_Vector     yprev;
     N_Vector     gy;
     N_Vector     fold;
     N_Vector     gold;
     N_Vector     delta;
     int          curiter;
     int          maxiters;
     long int     niters;
     long int     nconvfails;
     void        *ctest_data;
     int          print_level;
     FILE*        info_file;
   };

The following entries of the *content* field are always
allocated:

* ``Sys``        -- function for evaluating the nonlinear system,
* ``CTest``      -- function for checking convergence of the fixed point iteration,
* ``yprev``      -- ``N_Vector`` used to store previous fixed-point iterate,
* ``gy``         -- ``N_Vector`` used to store :math:`G(y)` in fixed-point algorithm,
* ``delta``      -- ``N_Vector`` used to store difference between successive fixed-point iterates,
* ``curiter``    -- the current number of iterations in the solve attempt,
* ``maxiters``   -- the maximum number of fixed-point iterations allowed in
  a solve,
* ``niters``     -- the total number of nonlinear iterations across all
  solves,
* ``nconvfails`` -- the total number of nonlinear convergence failures across all solves,
* ``ctest_data`` -- the data pointer passed to the convergence test function,
* ``m``          -- number of acceleration vectors,
* ``print_level`` - controls the amount of information to be printed to the info file, and
* ``info_file``   - the file where all informative (non-error) messages will be directed.

If Anderson acceleration is requested (i.e., :math:`m>0` in the call
to :c:func:`SUNNonlinSol_FixedPoint`), then the following items are also
allocated within the *content* field:

* ``imap``    -- index array used in acceleration algorithm (length ``m``),
* ``damping`` -- a flag indicating if damping is enabled,
* ``beta``    -- the damping parameter,
* ``R``       -- small matrix used in acceleration algorithm (length ``m*m``),
* ``gamma``   -- small vector used in acceleration algorithm (length ``m``),
* ``cvals``   -- small vector used in acceleration algorithm (length ``m+1``),
* ``df``      -- array of ``N_Vectors`` used in acceleration algorithm (length ``m``),
* ``dg``      -- array of ``N_Vectors`` used in acceleration algorithm (length ``m``),
* ``q``       -- array of ``N_Vectors`` used in acceleration algorithm (length ``m``),
* ``Xvecs``   -- ``N_Vector`` pointer array used in acceleration algorithm (length ``m+1``),
* ``fold``    -- ``N_Vector`` used in acceleration algorithm, and
* ``gold``    -- ``N_Vector`` used in acceleration algorithm.
