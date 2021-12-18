.. ----------------------------------------------------------------
   SUNDIALS Copyright Start
   Copyright (c) 2002-2021, Lawrence Livermore National Security
   and Southern Methodist University.
   All rights reserved.

   See the top-level LICENSE and NOTICE files for details.

   SPDX-License-Identifier: BSD-3-Clause
   SUNDIALS Copyright End
   ----------------------------------------------------------------

.. _IDAS.Mathematics:

***************************
Mathematical Considerations
***************************

IDAS solves the initial-value problem (IVP) for a DAE system of the general form

.. math::
   F(t,y,\dot{y}) = 0 \, , \quad y(t_0) = y_0 \, , \quad \dot{y}(t_0) =
   \dot{y}_0 \,
   :label: IDAS_DAE

where :math:`y`, :math:`\dot{y}`, and :math:`F` are vectors in :math:`{\bf
R}^N`, :math:`t` is the independent variable, :math:`\dot{y} = \mathrm
dy/\mathrm dt`, and initial values :math:`y_0`, :math:`\dot{y}_0` are given.
(Often :math:`t` is time, but it certainly need not be.)

Additionally, if :eq:`IDAS_DAE` depends on some parameters
:math:`p \in {\bf R}^{N_p}`, i.e.

.. math::
   \begin{split} & F(t, y, \dot y , p) = 0 \\ & y(t_0) = y_0(p) \, ,~ {\dot
   y}(t_0) = \dot y_0(p) \, , \end{split}
   :label: IDAS_DAE_p

IDAS can also compute first order derivative information, performing either
*forward sensitivity analysis* or *adjoint sensitivity analysis*. In the first
case, IDAS computes the sensitivities of the solution with respect to the
parameters :math:`p`, while in the second case, IDAS computes the gradient of a
*derived function* with respect to the parameters :math:`p`.

.. _IDAS.Mathematics.ivp_sol:

IVP solution
============

Prior to integrating a DAE initial-value problem, an important requirement is
that the pair of vectors :math:`y_0` and :math:`\dot{y}_0` are both initialized
to satisfy the DAE residual :math:`F(t_0,y_0, \dot{y}_0) = 0`.  For a class of
problems that includes so-called semi-explicit index-one systems, IDAS provides a
routine that computes consistent initial conditions from a user’s initial guess
:cite:p:`BHP:98`.  For this, the user must identify sub-vectors of :math:`y`
(not necessarily contiguous), denoted :math:`y_d` and :math:`y_a`, which are its
differential and algebraic parts, respectively, such that :math:`F` depends on
:math:`\dot{y}_d` but not on any components of :math:`\dot{y}_a`. The assumption
that the system is “index one” means that for a given :math:`t` and :math:`y_d`,
the system :math:`F(t,y,\dot{y}) = 0` defines :math:`y_a` uniquely. In this
case, a solver within IDAS computes :math:`y_a` and :math:`\dot{y}_d` at :math:`t
= t_0`, given :math:`y_d` and an initial guess for :math:`y_a`. A second
available option with this solver also computes all of :math:`y(t_0)` given
:math:`\dot{y}(t_0)`; this is intended mainly for quasi-steady-state problems,
where :math:`\dot{y}(t_0) = 0` is given.  In both cases, IDAS solves the system
:math:`F(t_0,y_0, \dot{y}_0) = 0` for the unknown components of :math:`y_0` and
:math:`\dot{y}_0`, using Newton iteration augmented with a line search global
strategy. In doing this, it makes use of the existing machinery that is to be
used for solving the linear systems during the integration, in combination with
certain tricks involving the step size (which is set artificially for this
calculation).  For problems that do not fall into either of these categories,
the user is responsible for passing consistent values, or risks failure in the
numerical integration.

The integration method used in IDAS is the variable-order, variable-coefficient
BDF (Backward Differentiation Formula), in fixed-leading-coefficient form
:cite:p:`BCP:96`.  The method order ranges from 1 to 5, with the BDF of order
:math:`q` given by the multistep formula

.. math::
   \sum_{i=0}^q \alpha_{n,i}y_{n-i} = h_n \dot{y}_n \, ,
   :label: IDAS_BDF

where :math:`y_n` and :math:`\dot{y}_n` are the computed approximations to
:math:`y(t_n)` and :math:`\dot{y}(t_n)`, respectively, and the step size is
:math:`h_n = t_n - t_{n-1}`.  The coefficients :math:`\alpha_{n,i}` are uniquely
determined by the order :math:`q`, and the history of the step sizes. The
application of the BDF :eq:`IDAS_BDF` to the DAE system :eq:`IDAS_DAE` results in a
nonlinear algebraic system to be solved at each step:

.. math::
   G(y_n) \equiv
   F \left( t_n , \, y_n , \,
      h_n^{-1} \sum_{i=0}^q \alpha_{n,i}y_{n-i} \right) = 0 \, .
   :label: IDAS_DAE_nls

By default IDAS solves :eq:`IDAS_DAE_nls` with a Newton iteration but IDAS also allows
for user-defined nonlinear solvers (see Chapter :numref:`SUNNonlinSol`). Each
Newton iteration requires the solution of a linear system of the form

.. math::
   J [y_{n(m+1)} - y_{n(m)}] = -G(y_{n(m)})  \, ,
   :label: IDAS_DAE_Newtoncorr

where :math:`y_{n(m)}` is the :math:`m`-th approximation to :math:`y_n`.  Here
:math:`J` is some approximation to the system Jacobian

.. math::
   J = \frac{\partial G}{\partial y}
   = \frac{\partial F}{\partial y} +
   \alpha\frac{\partial F}{\partial \dot{y}} \, ,
   :label: IDAS_DAE_Jacobian

where :math:`\alpha = \alpha_{n,0}/h_n`. The scalar :math:`\alpha` changes
whenever the step size or method order changes.

For the solution of the linear systems within the Newton iteration, IDAS provides
several choices, including the option of a user-supplied linear solver (see
Chapter :numref:`SUNLinSol`). The linear solvers distributed with SUNDIALS are
organized in two families, a *direct* family comprising direct linear solvers
for dense, banded, or sparse matrices and a *spils* family comprising scaled
preconditioned iterative (Krylov) linear solvers.  The methods offered through
these modules are as follows:

* dense direct solvers, including an internal implementation, an interface to
  BLAS/LAPACK, an interface to MAGMA :cite:p:`magma_ref` and an interface to
  the oneMKL library :cite:p:`oneAPI_site`,

* band direct solvers, including an internal implementation or an interface to BLAS/LAPACK,

* sparse direct solver interfaces to various libraries, including KLU
  :cite:p:`DaPa:10, KLU_site`, SuperLU_MT :cite:p:`Li:05,DGL:99,SuperLUMT_site`,
  SuperLU_Dist :cite:p:`GDL:07,LD:03,SLUUG:99,SuperLUDIST_site`, and cuSPARSE :cite:p:`cuSPARSE_site`,

* SPGMR, a scaled preconditioned GMRES (Generalized Minimal Residual method) solver,

* SPFGMR, a scaled preconditioned FGMRES (Flexible Generalized Minimal Residual method) solver,

* SPBCG, a scaled preconditioned Bi-CGStab (Bi-Conjugate Gradient Stable method) solver,

* SPTFQMR, a scaled preconditioned TFQMR (Transpose-Free Quasi-Minimal Residual method) solver, or

* PCG, a scaled preconditioned CG (Conjugate Gradient method) solver.

For large stiff systems, where direct methods are not feasible, the combination
of a BDF integrator and a preconditioned Krylov method yields a powerful tool
because it combines established methods for stiff integration, nonlinear
iteration, and Krylov (linear) iteration with a problem-specific treatment of
the dominant source of stiffness, in the form of the user-supplied
preconditioner matrix :cite:p:`BrHi:89`.  For the *spils* linear solvers with
IDAS, preconditioning is allowed only on the left (see
:numref:`IDAS.Mathematics.preconditioning`).  Note that the dense, band, and sparse
direct linear solvers can only be used with serial and threaded vector
representations.

In the process of controlling errors at various levels, IDAS uses a weighted
root-mean-square norm, denoted :math:`\|\cdot\|_{\mbox{WRMS}}`, for all
error-like quantities. The multiplicative weights used are based on the current
solution and on the relative and absolute tolerances input by the user, namely

.. math::
   W_i = \frac{1}{\text{rtol} \cdot |y_i| + \text{atol}_i }\, .
   :label: IDAS_errwt

Because :math:`1/W_i` represents a tolerance in the component :math:`y_i`, a
vector whose norm is 1 is regarded as “small.” For brevity, we will usually drop
the subscript WRMS on norms in what follows.

In the case of a matrix-based linear solver, the default Newton iteration is a
Modified Newton iteration, in that the Jacobian :math:`J` is fixed (and usually
out of date) throughout the nonlinear iterations, with a coefficient
:math:`\bar\alpha` in place of :math:`\alpha` in :math:`J`. However, in the case
that a matrix-free iterative linear solver is used, the default Newton iteration
is an Inexact Newton iteration, in which :math:`J` is applied in a matrix-free
manner, with matrix-vector products :math:`Jv` obtained by either difference
quotients or a user-supplied routine.  In this case, the linear residual
:math:`J\Delta y + G` is nonzero but controlled.  With the default Newton
iteration, the matrix :math:`J` and preconditioner matrix :math:`P` are updated
as infrequently as possible to balance the high costs of matrix operations
against other costs. Specifically, this matrix update occurs when:

* starting the problem,

* the value :math:`\bar\alpha` at the last update is such that :math:`\alpha /
  {\bar\alpha} < 3/5` or :math:`\alpha / {\bar\alpha} > 5/3`, or

* a non-fatal convergence failure occurred with an out-of-date :math:`J` or
  :math:`P`.

The above strategy balances the high cost of frequent matrix evaluations and
preprocessing with the slow convergence due to infrequent updates.  To reduce
storage costs on an update, Jacobian information is always reevaluated from
scratch.

The default stopping test for nonlinear solver iterations in IDAS ensures that
the iteration error :math:`y_n - y_{n(m)}` is small relative to :math:`y`
itself. For this, we estimate the linear convergence rate at all iterations
:math:`m>1` as

.. math:: R = \left( \frac{\delta_m}{\delta_1} \right)^{\frac{1}{m-1}} \, ,

where the :math:`\delta_m = y_{n(m)} - y_{n(m-1)}` is the correction at
iteration :math:`m=1,2,\ldots`. The nonlinear solver iteration is halted if
:math:`R>0.9`.  The convergence test at the :math:`m`-th iteration is then

.. math::
   S \| \delta_m \| < 0.33 \, ,
   :label: IDAS_DAE_nls_test

where :math:`S = R/(R-1)` whenever :math:`m>1` and :math:`R\le 0.9`. The user
has the option of changing the constant in the convergence test from its default
value of :math:`0.33`.  The quantity :math:`S` is set to :math:`S=20` initially
and whenever :math:`J` or :math:`P` is updated, and it is reset to :math:`S=100`
on a step with :math:`\alpha \neq \bar\alpha`.  Note that at :math:`m=1`, the
convergence test :eq:`IDAS_DAE_nls_test` uses an old value for :math:`S`. Therefore,
at the first nonlinear solver iteration, we make an additional test and stop the
iteration if :math:`\|\delta_1\| < 0.33 \cdot 10^{-4}` (since such a
:math:`\delta_1` is probably just noise and therefore not appropriate for use in
evaluating :math:`R`).  We allow only a small number (default value 4) of
nonlinear iterations.  If convergence fails with :math:`J` or :math:`P` current,
we are forced to reduce the step size :math:`h_n`, and we replace :math:`h_n` by
:math:`h_n/4`.  The integration is halted after a preset number (default
value 10) of convergence failures. Both the maximum number of allowable
nonlinear iterations and the maximum number of nonlinear convergence failures
can be changed by the user from their default values.

When an iterative method is used to solve the linear system, to minimize the
effect of linear iteration errors on the nonlinear and local integration error
controls, we require the preconditioned linear residual to be small relative to
the allowed error in the nonlinear iteration, i.e., :math:`\| P^{-1}(Jx+G) \| <
0.05 \cdot 0.33`.  The safety factor :math:`0.05` can be changed by the user.

When the Jacobian is stored using either the :ref:`SUNMATRIX_DENSE <SUNMatrix.Dense>`
or :ref:`SUNMATRIX_BAND <SUNMatrix.Band>`  matrix objects,
the Jacobian :math:`J` defined in :eq:`IDAS_DAE_Jacobian` can be either supplied by
the user or IDAS can compute :math:`J` internally by difference quotients. In the
latter case, we use the approximation

.. math::
   \begin{gathered}
     J_{ij} = [F_i(t,y+\sigma_j e_j,\dot{y}+\alpha\sigma_j e_j) -
               F_i(t,y,\dot{y})]/\sigma_j \, , \text{ with}\\
     \sigma_j = \sqrt{U} \max \left\{ |y_j|, |h\dot{y}_j|,1/W_j \right\}
                \mbox{sign}(h \dot{y}_j) \, ,\end{gathered}

where :math:`U` is the unit roundoff, :math:`h` is the current step size, and
:math:`W_j` is the error weight for the component :math:`y_j` defined by
:eq:`IDAS_errwt`.  We note that with sparse and user-supplied matrix objects,
the Jacobian *must* be supplied by a user routine.

In the case of an iterative linear solver, if a routine for :math:`Jv` is not
supplied, such products are approximated by

.. math:: Jv = [F(t,y+\sigma v,\dot{y}+\alpha\sigma v) - F(t,y,\dot{y})]/\sigma \, ,

where the increment :math:`\sigma = 1/\|v\|`. As an option, the user can specify
a constant factor that is inserted into this expression for :math:`\sigma`.

During the course of integrating the system, IDAS computes an estimate of the
local truncation error, LTE, at the :math:`n`-th time step, and requires this to
satisfy the inequality

.. math:: \| \mbox{LTE} \|_{\mbox{WRMS}} \leq 1 \, .

Asymptotically, LTE varies as :math:`h^{q+1}` at step size :math:`h` and order
:math:`q`, as does the predictor-corrector difference :math:`\Delta_n \equiv
y_n-y_{n(0)}`.  Thus there is a constant :math:`C` such that

.. math:: \mbox{LTE} = C \Delta_n + O(h^{q+2}) \, ,

and so the norm of LTE is estimated as :math:`|C| \cdot \|\Delta_n\|`.  In
addition, IDAS requires that the error in the associated polynomial interpolant
over the current step be bounded by 1 in norm. The leading term of the norm of
this error is bounded by :math:`\bar{C} \|\Delta_n\|` for another constant
:math:`\bar{C}`. Thus the local error test in IDAS is

.. math::
   \max\{ |C|, \bar{C} \} \|\Delta_n\| \leq 1 \, .
   :label: IDAS_lerrtest

A user option is available by which the algebraic components of the error vector
are omitted from the test :eq:`IDAS_lerrtest`, if these have been so identified.

In IDAS, the local error test is tightly coupled with the logic for selecting the
step size and order. First, there is an initial phase that is treated specially;
for the first few steps, the step size is doubled and the order raised (from its
initial value of 1) on every step, until (a) the local error test :eq:`IDAS_lerrtest`
fails, (b) the order is reduced (by the rules given below), or (c) the order
reaches 5 (the maximum). For step and order selection on the general step, IDAS
uses a different set of local error estimates, based on the asymptotic behavior
of the local error in the case of fixed step sizes.  At each of the orders
:math:`q'` equal to :math:`q`, :math:`q-1` (if :math:`q > 1`), :math:`q-2` (if
:math:`q > 2`), or :math:`q+1` (if :math:`q < 5`), there are constants
:math:`C(q')` such that the norm of the local truncation error at order
:math:`q'` satisfies

.. math:: \mbox{LTE}(q') = C(q') \| \phi(q'+1) \| + O(h^{q'+2}) \, ,

where :math:`\phi(k)` is a modified divided difference of order :math:`k` that
is retained by IDAS (and behaves asymptotically as :math:`h^k`).  Thus the local
truncation errors are estimated as ELTE\ :math:`(q') = C(q')\|\phi(q'+1)\|` to
select step sizes. But the choice of order in IDAS is based on the requirement
that the scaled derivative norms, :math:`\|h^k y^{(k)}\|`, are monotonically
decreasing with :math:`k`, for :math:`k` near :math:`q`. These norms are again
estimated using the :math:`\phi(k)`, and in fact

.. math:: \|h^{q'+1} y^{(q'+1)}\| \approx T(q') \equiv (q'+1) \mbox{ELTE}(q') \, .

The step/order selection begins with a test for monotonicity that is made even
*before* the local error test is performed. Namely, the order is reset to
:math:`q' = q-1` if (a) :math:`q=2` and :math:`T(1)\leq T(2)/2`, or (b) :math:`q
> 2` and :math:`\max\{T(q-1),T(q-2)\} \leq T(q)`; otherwise :math:`q' = q`. Next
the local error test :eq:`IDAS_lerrtest` is performed, and if it fails, the step is
redone at order :math:`q\leftarrow q'` and a new step size :math:`h'`. The
latter is based on the :math:`h^{q+1}` asymptotic behavior of
:math:`\mbox{ELTE}(q)`, and, with safety factors, is given by

.. math:: \eta = h'/h = 0.9/[2 \, \mbox{ELTE}(q)]^{1/(q+1)} \, .

The value of :math:`\eta` is adjusted so that :math:`0.25 \leq \eta \leq 0.9`
before setting :math:`h \leftarrow h' = \eta h`. If the local error test fails a
second time, IDAS uses :math:`\eta = 0.25`, and on the third and subsequent
failures it uses :math:`q = 1` and :math:`\eta = 0.25`. After 10 failures, IDAS
returns with a give-up message.

As soon as the local error test has passed, the step and order for the next step
may be adjusted. No such change is made if :math:`q' = q-1` from the prior test,
if :math:`q = 5`, or if :math:`q` was increased on the previous step. Otherwise,
if the last :math:`q+1` steps were taken at a constant order :math:`q < 5` and a
constant step size, IDAS considers raising the order to :math:`q+1`. The logic is
as follows: (a) If :math:`q = 1`, then reset :math:`q = 2` if :math:`T(2) <
T(1)/2`. (b) If :math:`q > 1` then

* reset :math:`q \leftarrow q-1` if :math:`T(q-1) \leq \min\{T(q),T(q+1)\}`;

* else reset :math:`q \leftarrow q+1` if :math:`T(q+1) < T(q)`;

* leave :math:`q` unchanged otherwise :math:`[`\ then :math:`T(q-1) > T(q) \leq
  T(q+1)]`.

In any case, the new step size :math:`h'` is set much as before:

.. math:: \eta = h'/h = 1/[2 \, \mbox{ELTE}(q)]^{1/(q+1)} \, .

The value of :math:`\eta` is adjusted such that (a) if :math:`\eta > 2`,
:math:`\eta` is reset to 2; (b) if :math:`\eta \leq 1`, :math:`\eta` is
restricted to :math:`0.5 \leq \eta \leq 0.9`; and (c) if :math:`1 < \eta < 2` we
use :math:`\eta = 1`.  Finally :math:`h` is reset to :math:`h' = \eta h`. Thus
we do not increase the step size unless it can be doubled. See :cite:p:`BCP:96`
for details.

IDAS permits the user to impose optional inequality constraints on individual
components of the solution vector :math:`y`. Any of the following four
constraints can be imposed: :math:`y_i > 0`, :math:`y_i < 0`, :math:`y_i \geq
0`, or :math:`y_i \leq 0`.  The constraint satisfaction is tested after a
successful nonlinear system solution.  If any constraint fails, we declare a
convergence failure of the nonlinear iteration and reduce the step size. Rather
than cutting the step size by some arbitrary factor, IDAS estimates a new step
size :math:`h'` using a linear approximation of the components in :math:`y` that
failed the constraint test (including a safety factor of :math:`0.9` to cover
the strict inequality case). These additional constraints are also imposed
during the calculation of consistent initial conditions.  If a step fails to
satisfy the constraints repeatedly within a step attempt then the integration is
halted and an error is returned. In this case the user may need to employ other
strategies as discussed in :numref:`IDAS.Usage.SIM.user_callable.idatolerances` to
satisfy the inequality constraints.

Normally, IDAS takes steps until a user-defined output value :math:`t =
t_{\mbox{out}}` is overtaken, and then computes :math:`y(t_{\mbox{out}})` by
interpolation. However, a “one step” mode option is available, where control
returns to the calling program after each step. There are also options to force
IDAS not to integrate past a given stopping point :math:`t = t_{\mbox{stop}}`.

.. _IDAS.Mathematics.preconditioning:

Preconditioning
===============

When using a nonlinear solver that requires the solution of a linear system of
the form :math:`J \Delta y = - G` (e.g., the default Newton iteration), IDAS
makes repeated use of a linear solver.  If this linear system solve is done with
one of the scaled preconditioned iterative linear solvers supplied with
SUNDIALS, these solvers are rarely successful if used without preconditioning;
it is generally necessary to precondition the system in order to obtain
acceptable efficiency.  A system :math:`A x = b` can be preconditioned on the
left, on the right, or on both sides. The Krylov method is then applied to a
system with the matrix :math:`P^{-1}A`, or :math:`AP^{-1}`, or :math:`P_L^{-1} A
P_R^{-1}`, instead of :math:`A`.  However, within IDAS, preconditioning is
allowed *only* on the left, so that the iterative method is applied to systems
:math:`(P^{-1}J)\Delta y = -P^{-1}G`.  Left preconditioning is required to make
the norm of the linear residual in the nonlinear iteration meaningful; in
general, :math:`\| J \Delta y + G \|` is meaningless, since the weights used in
the WRMS-norm correspond to :math:`y`.

In order to improve the convergence of the Krylov iteration, the preconditioner
matrix :math:`P` should in some sense approximate the system matrix :math:`A`.
Yet at the same time, in order to be cost-effective, the matrix :math:`P` should
be reasonably efficient to evaluate and solve. Finding a good point in this
tradeoff between rapid convergence and low cost can be very difficult. Good
choices are often problem-dependent (for example, see :cite:p:`BrHi:89` for an
extensive study of preconditioners for reaction-transport systems).

Typical preconditioners used with IDAS are based on approximations to the
iteration matrix of the systems involved; in other words, :math:`P \approx
{\partial F}/{\partial y} + \alpha\, {\partial F}/{\partial \dot{y}}`,
where :math:`\alpha` is a scalar inversely proportional to the integration step
size :math:`h`.  Because the Krylov iteration occurs within a nonlinear solver
iteration and further also within a time integration, and since each of these
iterations has its own test for convergence, the preconditioner may use a very
crude approximation, as long as it captures the dominant numerical feature(s) of
the system. We have found that the combination of a preconditioner with the
Newton-Krylov iteration, using even a fairly poor approximation to the Jacobian,
can be surprisingly superior to using the same matrix without Krylov
acceleration (i.e., a modified Newton iteration), as well as to using the
Newton-Krylov method with no preconditioning.

.. _IDAS.Mathematics.rootfinding:

Rootfinding
===========

The IDAS solver has been augmented to include a rootfinding feature. This means
that, while integratnuming the Initial Value Problem :eq:`IDAS_DAE`, IDAS can also
find the roots of a set of user-defined functions :math:`g_i(t,y,\dot{y})` that
depend on :math:`t`, the solution vector :math:`y = y(t)`, and its :math:`t-`\
derivative :math:`\dot{y}(t)`. The number of these root functions is arbitrary,
and if more than one :math:`g_i` is found to have a root in any given interval,
the various root locations are found and reported in the order that they occur
on the :math:`t` axis, in the direction of integration.

Generally, this rootfinding feature finds only roots of odd multiplicity,
corresponding to changes in sign of :math:`g_i(t,y(t),\dot{y}(t))`, denoted
:math:`g_i(t)` for short. If a user root function has a root of even
multiplicity (no sign change), it will probably be missed by IDAS. If such a root
is desired, the user should reformulate the root function so that it changes
sign at the desired root.

The basic scheme used is to check for sign changes of any :math:`g_i(t)` over
each time step taken, and then (when a sign change is found) to home in on the
root (or roots) with a modified secant method :cite:p:`HeSh:80`.  In addition,
each time :math:`g` is computed, IDAS checks to see if :math:`g_i(t) = 0`
exactly, and if so it reports this as a root. However, if an exact zero of any
:math:`g_i` is found at a point :math:`t`, IDAS computes :math:`g` at :math:`t +
\delta` for a small increment :math:`\delta`, slightly further in the direction
of integration, and if any :math:`g_i(t + \delta)=0` also, IDAS stops and reports
an error. This way, each time IDAS takes a time step, it is guaranteed that the
values of all :math:`g_i` are nonzero at some past value of :math:`t`, beyond
which a search for roots is to be done.

At any given time in the course of the time-stepping, after suitable checking
and adjusting has been done, IDAS has an interval :math:`(t_{lo},t_{hi}]` in
which roots of the :math:`g_i(t)` are to be sought, such that :math:`t_{hi}` is
further ahead in the direction of integration, and all :math:`g_i(t_{lo}) \neq
0`. The endpoint :math:`t_{hi}` is either :math:`t_n`, the end of the time step
last taken, or the next requested output time :math:`t_{\mbox{out}}` if this
comes sooner. The endpoint :math:`t_{lo}` is either :math:`t_{n-1}`, or the last
output time :math:`t_{\mbox{out}}` (if this occurred within the last step), or
the last root location (if a root was just located within this step), possibly
adjusted slightly toward :math:`t_n` if an exact zero was found. The algorithm
checks :math:`g` at :math:`t_{hi}` for zeros and for sign changes in
:math:`(t_{lo},t_{hi})`. If no sign changes are found, then either a root is
reported (if some :math:`g_i(t_{hi}) = 0`) or we proceed to the next time
interval (starting at :math:`t_{hi}`). If one or more sign changes were found,
then a loop is entered to locate the root to within a rather tight tolerance,
given by

.. math:: \tau = 100 * U * (|t_n| + |h|)~~~ (U = \mbox{unit roundoff}) ~.

Whenever sign changes are seen in two or more root functions, the one deemed
most likely to have its root occur first is the one with the largest value of
:math:`|g_i(t_{hi})|/|g_i(t_{hi}) - g_i(t_{lo})|`, corresponding to the closest
to :math:`t_{lo}` of the secant method values.  At each pass through the loop, a
new value :math:`t_{mid}` is set, strictly within the search interval, and the
values of :math:`g_i(t_{mid})` are checked. Then either :math:`t_{lo}` or
:math:`t_{hi}` is reset to :math:`t_{mid}` according to which subinterval is
found to have the sign change. If there is none in :math:`(t_{lo},t_{mid})` but
some :math:`g_i(t_{mid}) = 0`, then that root is reported. The loop continues
until :math:`|t_{hi}-t_{lo}| < \tau`, and then the reported root location is
:math:`t_{hi}`.

In the loop to locate the root of :math:`g_i(t)`, the formula for
:math:`t_{mid}` is

.. math::

   t_{mid} = t_{hi} - (t_{hi} - t_{lo})
                g_i(t_{hi}) / [g_i(t_{hi}) - \alpha g_i(t_{lo})] ~,

where :math:`\alpha` a weight parameter. On the first two passes through the
loop, :math:`\alpha` is set to :math:`1`, making :math:`t_{mid}` the secant
method value. Thereafter, :math:`\alpha` is reset according to the side of the
subinterval (low vs high, i.e. toward :math:`t_{lo}` vs toward :math:`t_{hi}`)
in which the sign change was found in the previous two passes. If the two sides
were opposite, :math:`\alpha` is set to 1. If the two sides were the same,
:math:`\alpha` is halved (if on the low side) or doubled (if on the high
side). The value of :math:`t_{mid}` is closer to :math:`t_{lo}` when
:math:`\alpha < 1` and closer to :math:`t_{hi}` when :math:`\alpha > 1`. If the
above value of :math:`t_{mid}` is within :math:`\tau/2` of :math:`t_{lo}` or
:math:`t_{hi}`, it is adjusted inward, such that its fractional distance from
the endpoint (relative to the interval size) is between 0.1 and 0.5 (0.5 being
the midpoint), and the actual distance from the endpoint is at least
:math:`\tau/2`.

.. _IDAS.Mathematics.Purequad:

Pure quadrature integration
===========================

In many applications, and most notably during the backward integration phase of
an adjoint sensitivity analysis run :numref:`IDAS.Mathematics.ASA` it is of
interest to compute integral quantities of the form

.. math::
    z(t) = \int_{t_0}^t q(\tau, y(\tau), \dot{y}(\tau), p) \, d\tau \, .
    :label: IDAS_QUAD

The most effective approach to compute :math:`z(t)` is to extend the original
problem with the additional ODEs (obtained by applying Leibnitz’s
differentiation rule):

.. math:: \dot z = q(t,y,\dot{y},p) \, , \quad z(t_0) = 0 \, .

Note that this is equivalent to using a quadrature method based on the
underlying linear multistep polynomial representation for :math:`y(t)`.

This can be done at the “user level” by simply exposing to IDAS the extended DAE
system :eq:`IDAS_DAE_p` + :eq:`IDAS_QUAD`. However, in the context of an implicit
integration solver, this approach is not desirable since the nonlinear solver
module will require the Jacobian (or Jacobian-vector product) of this extended
DAE. Moreover, since the additional states, :math:`z`, do not enter the
right-hand side of the ODE :eq:`IDAS_QUAD` and therefore the residual of the extended
DAE system does not depend on :math:`z`, it is much more efficient to treat the
ODE system :eq:`IDAS_QUAD` separately from the original DAE system :eq:`IDAS_DAE_p` by
“taking out” the additional states :math:`z` from the nonlinear system
:eq:`IDAS_DAE_nls` that must be solved in the correction step of the LMM. Instead,
“corrected” values :math:`z_n` are computed explicitly as

.. math::
   z_n = \frac{1}{\alpha_{n,0}} \left(
       h_n q(t_n, y_n, \dot{y}_n, p) - \sum_{i=1}^q \alpha_{n,i} z_{n-i} \right)
       \, ,

once the new approximation :math:`y_n` is available.

The quadrature variables :math:`z` can be optionally included in the error test,
in which case corresponding relative and absolute tolerances must be provided.

.. _IDAS.Mathematics.FSA:

Forward sensitivity analysis
============================

Typically, the governing equations of complex, large-scale models depend on
various parameters, through the right-hand side vector and/or through the vector
of initial conditions, as in :eq:`IDAS_DAE_p`. In addition to numerically solving the
DAEs, it may be desirable to determine the sensitivity of the results with
respect to the model parameters. Such sensitivity information can be used to
estimate which parameters are most influential in affecting the behavior of the
simulation or to evaluate optimization gradients (in the setting of dynamic
optimization, parameter estimation, optimal control, etc.).

The *solution sensitivity* with respect to the model parameter :math:`p_i` is
defined as the vector :math:`s_i (t) = {\partial y(t)}/{\partial p_i}` and
satisfies the following *forward sensitivity equations* (or *sensitivity
equations* for short):

.. math::
   \begin{split}
   & \frac{\partial F}{\partial y} s_i + \frac{\partial F}{\partial \dot y} {\dot s_i} + \frac{\partial F}{\partial p_i} = 0\\
   & s_i(t_0) = \frac{\partial y_{0}(p)}{\partial p_i} \, ,~ \dot s_i(t_0) =  \frac{\partial \dot y_{0}(p)}{\partial p_i} \, ,
   \end{split}
  :label: IDAS_sens_eqns

obtained by applying the chain rule of differentiation to the original
DAEs :eq:`IDAS_DAE_p`.

When performing forward sensitivity analysis, IDAS carries out the time
integration of the combined system, :eq:`IDAS_DAE_p` and :eq:`IDAS_sens_eqns`, by viewing
it as a DAE system of size :math:`N(N_s+1)`, where :math:`N_s` is the number of
model parameters :math:`p_i`, with respect to which sensitivities are desired
(:math:`N_s \le N_p`). However, major improvements in efficiency can be made by
taking advantage of the special form of the sensitivity equations as
linearizations of the original DAEs. In particular, the original DAE system and
all sensitivity systems share the same Jacobian matrix :math:`J` in
:eq:`IDAS_DAE_Jacobian`.

The sensitivity equations are solved with the same linear multistep formula that
was selected for the original DAEs and the same linear solver is used in the
correction phase for both state and sensitivity variables. In addition, IDAS
offers the option of including (*full error control*) or excluding (*partial
error control*) the sensitivity variables from the local error test.

Forward sensitivity methods
---------------------------

In what follows we briefly describe three methods that have been proposed for the solution of the combined DAE and
sensitivity system for the vector :math:`{\hat y} = [y, s_1, \ldots , s_{N_s}]`.

-  *Staggered Direct* In this approach :cite:p:`CaSt:85`, the nonlinear system
   :eq:`IDAS_DAE_nls` is first solved and, once an acceptable numerical solution is
   obtained, the sensitivity variables at the new step are found by directly
   solving :eq:`IDAS_sens_eqns` after the BDF discretization is used to eliminate
   :math:`{\dot s}_i`. Although the system matrix of the above linear system is
   based on exactly the same information as the matrix :math:`J` in
   :eq:`IDAS_DAE_Jacobian`, it must be updated and factored at every step of the
   integration, in contrast to an evaluation of :math:`J` which is updated only
   occasionally. For problems with many parameters (relative to the problem
   size), the staggered direct method can outperform the methods described
   below :cite:p:`LPZ:99`. However, the computational cost associated with matrix
   updates and factorizations makes this method unattractive for problems with
   many more states than parameters (such as those arising from
   semidiscretization of PDEs) and is therefore not implemented in IDAS.

-  *Simultaneous Corrector* In this method :cite:p:`MaPe:97`, the discretization
   is applied simultaneously to both the original equations :eq:`IDAS_DAE_p` and the
   sensitivity systems :eq:`IDAS_sens_eqns` resulting in an “extended” nonlinear
   system :math:`{\hat G}({\hat y}_n) = 0` where :math:`{\hat y_n} = [ y_n,
   \ldots, s_i, \ldots ]`. This combined nonlinear system can be solved using a
   modified Newton method as in :eq:`IDAS_DAE_Newtoncorr` by solving the corrector
   equation

   .. math::
      {\hat J}[{\hat y}_{n(m+1)}-{\hat y}_{n(m)}]=-{\hat G}({\hat y}_{n(m)})
      :label: IDAS_Newton_sim

   at each iteration, where

   .. math::
      {\hat J} =
          \begin{bmatrix}
            J       &        &        &        &   \\
            J_1     & J      &        &        &   \\
            J_2     & 0      & J      &        &   \\
            \vdots  & \vdots & \ddots & \ddots &   \\
            J_{N_s} & 0      & \ldots & 0      & J
          \end{bmatrix} \, ,

   :math:`J` is defined as in :eq:`IDAS_DAE_Jacobian`, and :math:`J_i =
   ({\partial}/{\partial y})\left[ F_y s_i + F_{\dot y} {\dot s_i} + F_{p_i}
   \right]`. It can be shown that 2-step quadratic convergence can be retained
   by using only the block-diagonal portion of :math:`{\hat J}` in the corrector
   equation :eq:`IDAS_Newton_sim`. This results in a decoupling that allows the reuse
   of :math:`J` without additional matrix factorizations. However, the sum
   :math:`F_y s_i + F_{\dot y} {\dot s_i} + F_{p_i}` must still be reevaluated
   at each step of the iterative process :eq:`IDAS_Newton_sim` to update the
   sensitivity portions of the residual :math:`{\hat G}`.

-  *Staggered corrector* In this approach :cite:p:`FTB:97`, as in the staggered
   direct method, the nonlinear system :eq:`IDAS_DAE_nls` is solved first using the
   Newton iteration :eq:`IDAS_DAE_Newtoncorr`. Then, for each sensitivity vector
   :math:`\xi \equiv s_i`, a separate Newton iteration is used to solve the
   sensitivity system :eq:`IDAS_sens_eqns`:

   .. math::
      \begin{gathered}
          J [\xi_{n(m+1)} - \xi_{n(m)}]= \\
          - \left[
            F_y (t_n, y_n, \dot y_n) \xi_{n(m)}
            + F_{\dot y} (t_n, y_n, \dot y_n) \cdot
            h_n^{-1} \left(
              \alpha_{n,0} \xi_{n(m)} + \sum_{i=1}^q \alpha_{n,i} \xi_{n-i}
            \right)
            + F_{p_i} (t_n, y_n, \dot y_n)
          \right] \, .
      \end{gathered}
      :label: IDAS_stgr_iterations

   In other words, a modified Newton iteration is used to solve a linear system.
   In this approach, the matrices :math:`\partial F / \partial y`, :math:`\partial F / \partial \dot y` and vectors
   :math:`\partial f / \partial p_i` need be updated only once per integration step, after the
   state correction phase :eq:`IDAS_DAE_Newtoncorr` has converged.

IDAS implements both the simultaneous corrector method and the staggered
corrector method.

An important observation is that the staggered corrector method, combined with a
Krylov linear solver, effectively results in a staggered direct method. Indeed,
the Krylov solver requires only the action of the matrix :math:`J` on a vector,
and this can be provided with the current Jacobian information. Therefore, the
modified Newton procedure :eq:`IDAS_stgr_iterations` will theoretically converge
after one iteration.

Selection of the absolute tolerances for sensitivity variables
--------------------------------------------------------------

If the sensitivities are included in the error test, IDAS provides an automated
estimation of absolute tolerances for the sensitivity variables based on the
absolute tolerance for the corresponding state variable. The relative tolerance
for sensitivity variables is set to be the same as for the state variables. The
selection of absolute tolerances for the sensitivity variables is based on the
observation that the sensitivity vector :math:`s_i` will have units of
:math:`[y]/[p_i]`. With this, the absolute tolerance for the :math:`j`-th
component of the sensitivity vector :math:`s_i` is set to
:math:`{\mbox{atol}_j}/{|{\bar p}_i|}`, where :math:`\mbox{atol}_j`
are the absolute tolerances for the state variables and :math:`\bar p` is a
vector of scaling factors that are dimensionally consistent with the model
parameters :math:`p` and give an indication of their order of magnitude. This
choice of relative and absolute tolerances is equivalent to requiring that the
weighted root-mean-square norm of the sensitivity vector :math:`s_i` with
weights based on :math:`s_i` be the same as the weighted root-mean-square norm
of the vector of scaled sensitivities :math:`{\bar s}_i = |{\bar p}_i| s_i` with
weights based on the state variables (the scaled sensitivities :math:`{\bar
s}_i` being dimensionally consistent with the state variables). However, this
choice of tolerances for the :math:`s_i` may be a poor one, and the user of IDAS
can provide different values as an option.

Evaluation of the sensitivity right-hand side
---------------------------------------------

There are several methods for evaluating the residual functions in the
sensitivity systems :eq:`IDAS_sens_eqns`: analytic evaluation, automatic
differentiation, complex-step approximation, and finite differences (or
directional derivatives). IDAS provides all the software hooks for implementing
interfaces to automatic differentiation (AD) or complex-step approximation;
future versions will include a generic interface to AD-generated functions. At
the present time, besides the option for analytical sensitivity right-hand sides
(user-provided), IDAS can evaluate these quantities using various finite
difference-based approximations to evaluate the terms
:math:`(\partial F / \partial y) s_i + (\partial F / \partial \dot y) \dot s_i`
and :math:`(\partial f / \partial p_i)`, or using directional derivatives to
evaluate :math:`\left[ (\partial F / \partial y) s_i + (\partial F / \partial \dot y) \dot s_i + (\partial f / \partial p_i) \right]`.
As is typical for finite differences, the proper choice of perturbations is a
delicate matter. IDAS takes into account several problem-related features: the
relative DAE error tolerance :math:`\mbox{rtol}`, the machine unit roundoff
:math:`U`, the scale factor :math:`{\bar p}_i`, and the weighted
root-mean-square norm of the sensitivity vector :math:`s_i`.

Using central finite differences as an example, the two terms
:math:`(\partial F / \partial y) s_i + (\partial F / \partial \dot y) \dot s_i` and :math:`\partial f / \partial p_i`
in :eq:`IDAS_sens_eqns` can be evaluated either separately:

.. math::
   \frac{\partial F}{\partial y} s_i + \frac{\partial F}{\partial y}p \dot s_i \approx
   \frac{F(t, y+\sigma_y s_i, \dot y + \sigma_y \dot s_i , p)- F(t, y-\sigma_y s_i, \dot y - \sigma_y \dot s_i , p)}{2\,\sigma_y} \, ,
   :label: IDAS_fd2

.. math::
   \frac{\partial F}{\partial p_i} \approx
   \frac{F(t, y, \dot y, p + \sigma_i e_i)- F(t, y, \dot y, p - \sigma_i e_i)}{2\,\sigma_i} \, ,
   :label: IDAS_fd2p

.. math::
   \sigma_i = |{\bar p}_i| \sqrt{\max( \mbox{rtol} , U)} \, , \quad
   \sigma_y = \frac{1}{\max(1/\sigma_i, \|s_i\|_{\mbox{WRMS}}/|{\bar p}_i|)} \,

or simultaneously:

.. math::
   \frac{\partial F}{\partial y} s_i + \frac{\partial F}{\partial y}p \dot s_i + \frac{\partial F}{\partial p_i} \approx
   \frac{F(t, y+\sigma s_i, \dot y + \sigma \dot s_i , p + \sigma e_i) -
   F(t, y-\sigma s_i, \dot y - \sigma \dot s_i , p - \sigma e_i)}{2\,\sigma} \, ,
   :label: IDAS_dd2

.. math::
   \sigma = \min(\sigma_i, \sigma_y) \, ,

or by adaptively switching between :eq:`IDAS_fd2` + :eq:`IDAS_fd2p` and :eq:`IDAS_dd2`,
depending on the relative size of the two finite difference increments
:math:`\sigma_i` and :math:`\sigma_y`. In the adaptive scheme, if :math:`\rho =
\max(\sigma_i/\sigma_y,\sigma_y/\sigma_i)`, we use separate evaluations if
:math:`\rho > \rho_{max}` (an input value), and simultaneous evaluations otherwise.

These procedures for choosing the perturbations (:math:`\sigma_i`,
:math:`\sigma_y`, :math:`\sigma`) and switching between derivative formulas have
also been implemented for one-sided difference formulas. Forward finite
differences can be applied to :math:`(\partial F / \partial y) s_i + (\partial F / \partial \dot y) \dot s_i` and
:math:`\partial F / \partial p_i` separately, or the single directional derivative formula

.. math::

   \frac{\partial F}{\partial y} s_i + \frac{\partial F}{\partial y}p \dot s_i + \frac{\partial F}{\partial p_i} \approx
   \frac{F(t, y+\sigma s_i, \dot y + \sigma \dot s_i , p + \sigma e_i) - F(t, y, \dot y, p)}{\sigma}

can be used. In IDAS, the default value of :math:`\rho_{max}=0` indicates the use
of the second-order centered directional derivative formula :eq:`IDAS_dd2`
exclusively. Otherwise, the magnitude of :math:`\rho_{max}` and its sign (positive
or negative) indicates whether this switching is done with regard to (centered
or forward) finite differences, respectively.

Quadratures depending on forward sensitivities
----------------------------------------------

If pure quadrature variables are also included in the problem definition (see
:numref:`IDAS.Mathematics.Purequad`), IDAS does *not* carry their sensitivities
automatically. Instead, we provide a more general feature through which
integrals depending on both the states :math:`y` of :eq:`IDAS_DAE_p` and the state
sensitivities :math:`s_i` of :eq:`IDAS_sens_eqns` can be evaluated. In other words,
IDAS provides support for computing integrals of the form:

.. math::

   \bar z(t) = \int_{t_0}^t \bar q(\tau, y(\tau), \dot{y}(\tau), s_1(\tau), \ldots,
                 s_{N_p}(\tau),p) \, \mathrm d\tau \, .

If the sensitivities of the quadrature variables :math:`z` of :eq:`IDAS_QUAD` are
desired, these can then be computed by using:

.. math:: \bar q_i = q_y s_i + q_{\dot{y}} \dot{s}_i + q_{p_i} \, , \quad i = 1,\ldots,N_p \, ,

as integrands for :math:`\bar z`, where :math:`q_y`, :math:`q_{\dot{y}}`, and
:math:`q_p` are the partial derivatives of the integrand function :math:`q` of
:eq:`IDAS_QUAD`.

As with the quadrature variables :math:`z`, the new variables :math:`\bar z` are
also excluded from any nonlinear solver phase and “corrected” values :math:`\bar
z_n` are obtained through explicit formulas.


.. _IDAS.Mathematics.ASA:

Adjoint sensitivity analysis
============================

In the *forward sensitivity approach* described in the previous section,
obtaining sensitivities with respect to :math:`N_s` parameters is roughly
equivalent to solving an DAE system of size :math:`(1+N_s) N`. This can become
prohibitively expensive, especially for large-scale problems, if sensitivities
with respect to many parameters are desired. In this situation, the *adjoint
sensitivity method* is a very attractive alternative, provided that we do not
need the solution sensitivities :math:`s_i`, but rather the gradients with
respect to model parameters of a relatively few derived functionals of the
solution. In other words, if :math:`y(t)` is the solution of :eq:`IDAS_DAE_p`, we
wish to evaluate the gradient :math:`{\mathrm dG}/{\mathrm dp}` of

.. math::
   G(p) = \int_{t_0}^T g(t, y, p) \mathrm dt \, ,
   :label: IDAS_G

or, alternatively, the gradient :math:`{\mathrm dg}/{\mathrm dp}` of the function
:math:`g(t, y,p)` at the final time :math:`t = T`. The function :math:`g` must be
smooth enough that :math:`\partial g / \partial y` and :math:`\partial g / \partial p`
exist and are bounded.

In what follows, we only sketch the analysis for the sensitivity problem for
both :math:`G` and :math:`g`. For details on the derivation see
:cite:p:`CLPS:03`.

Sensitivity of :math:`G(p)`
---------------------------

We focus first on solving the sensitivity problem for :math:`G(p)` defined by :eq:`IDAS_G`. Introducing a Lagrange
multiplier :math:`\lambda`, we form the augmented objective function

.. math:: I(p) = G(p) - \int_{t_0}^T \lambda^*F(t, y,\dot y, p) \mathrm dt.

Since :math:`F(t, y,\dot y, p)=0`, the sensitivity of :math:`G` with respect to :math:`p` is

.. math::
   \frac{dG}{dp} = \frac{dI}{dp}
   =\int_{t_0}^T(g_p + g_yy_p)\mathrm dt - \int_{t_0}^T \lambda^*( F_p + F_yy_p +
   F_{\dot{y}}\dot{y}_p)\mathrm dt,
   :label: IDAS_dGdp1

where subscripts on functions such as :math:`F` or :math:`g` are used to denote
partial derivatives. By integration by parts, we have

.. math::
   \int_{t_0}^T \lambda^* F_{\dot{y}} \dot{y}_p \mathrm dt =
     (\lambda^* F_{\dot{y}}y_p) \big\vert_{t_0}^{T}
     - \int_{t_0}^T (\lambda^* F_{\dot{y}})' y_p \mathrm dt ,

where :math:`(\cdots)'` denotes the :math:`t-`\ derivative. Thus equation
:eq:`IDAS_dGdp1` becomes

.. math::
   \frac{dG}{dp} = \int_{t_0}^T \left(g_p - \lambda^*F_p \right) \mathrm dt -
       \int_{t_0}^T \left[-g_y + \lambda^*F_y - (\lambda^*F_{\dot y})'\right]y_p \mathrm dt
        - (\lambda^* F_{\dot{y}} y_p) \big\vert_{t_0}^{T}.

Now by requiring :math:`\lambda` to satisfy

.. math::
   (\lambda^*F_{\dot{y}})' - \lambda^*F_y = -g_y ,
   :label: IDAS_adj_eqns

we obtain

.. math::
   \frac{dG}{dp} = \int_{t_0}^T \left(
     g_p - \lambda^*F_p \right) \mathrm dt - (\lambda^* F_{\dot{y}}y_p)\big\vert_{t_0}^T .
   :label: IDAS_dGdp

Note that :math:`y_p` at :math:`t=t_0` is the sensitivity of the initial
conditions with respect to :math:`p`, which is easily obtained. To find the
initial conditions (at :math:`t = T`) for the adjoint system, we must take into
consideration the structure of the DAE system.

For index-0 and index-1 DAE systems, we can simply take

.. math::
   \lambda^*F_{\dot y}\big\vert_{t=T} = 0,
   :label: IDAS_ad-init1

yielding the sensitivity equation for :math:`{dG}/{dp}`

.. math::
   \frac{dG}{dp} = \int_{t_0}^T \left(
     g_p - \lambda^*F_p \right) \mathrm dt
   + (\lambda^* F_{\dot{y}}y_p)\big\vert_{t=t_0} .
   :label: IDAS_sensi12

This choice will not suffice for a Hessenberg index-2 DAE system. For a
derivation of proper final conditions in such cases, see :cite:p:`CLPS:03`.

The first thing to notice about the adjoint system :eq:`IDAS_adj_eqns` is that there
is no explicit specification of the parameters :math:`p`; this implies that,
once the solution :math:`\lambda` is found, the formula :eq:`IDAS_dGdp` can then be
used to find the gradient of :math:`G` with respect to any of the parameters
:math:`p`. The second important remark is that the adjoint system :eq:`IDAS_adj_eqns`
is a terminal value problem which depends on the solution :math:`y(t)` of the
original IVP :eq:`IDAS_DAE_p`. Therefore, a procedure is needed for providing the
states :math:`y` obtained during a forward integration phase of :eq:`IDAS_DAE_p` to
IDAS during the backward integration phase of :eq:`IDAS_adj_eqns`. The approach
adopted in IDAS, based on *checkpointing*, is described in
:numref:`IDAS.Mathematics.ASA.Checkpointing` below.

Sensitivity of :math:`g(T,p)`
-----------------------------

Now let us consider the computation of :math:`{\mathrm dg}/{\mathrm dp}(T)`. From
:math:`{\mathrm dg}/{\mathrm dp}(T) = ({\mathrm d}/{\mathrm dT})({\mathrm dG}/{\mathrm dp})` and
equation :eq:`IDAS_dGdp`, we have

.. math::
   \frac{\mathrm dg}{\mathrm dp} = (g_p - \lambda^*F_p)(T) - \int_{t_0}^T \lambda^*_TF_p \mathrm dt +
    (\lambda^*_T F_{\dot{y}}y_p)\bigg\vert_{t=t_0} - \frac{\mathrm d(\lambda^*F_{\dot y}y_p)}{\mathrm dT}
   :label: IDAS_dlowgdp1

where :math:`\lambda_T` denotes :math:`{\partial \lambda}/{\partial T}`. For
index-0 and index-1 DAEs, we obtain

.. math:: \frac{\mathrm d(\lambda^*F_{\dot y}y_p)\big\vert_{t=T}}{\mathrm dT} = 0 ,

while for a Hessenberg index-2 DAE system we have

.. math::
   \frac{\mathrm d(\lambda^*F_{\dot y}y_p)\big\vert_{t=T}}{\mathrm dT} =
   -\left.\frac{\mathrm d(g_{y^a}(CB)^{-1}f^2_p)}{\mathrm dt}\right|_{t=T} .

The corresponding adjoint equations are

.. math::
   (\lambda^*_TF_{\dot y})'  - \lambda^*_T F_y = 0.
   :label: IDAS_adj1_eqns

For index-0 and index-1 DAEs (as shown above, the index-2 case is different), to
find the boundary condition for this equation we write :math:`\lambda` as
:math:`\lambda(t, T)` because it depends on both :math:`t` and :math:`T`. Then

.. math:: \lambda^*(T, T) F_{\dot{y}}\bigg\vert_{t=T}  = 0.

Taking the total derivative, we obtain

.. math::
   (\lambda_t + \lambda_T)^*(T, T) F_{\dot{y}}\bigg\vert_{t=T}  +
   \lambda^*(T,T)\frac{\mathrm dF_{\dot{y}}}{\mathrm dt}\bigg\vert_{t=T} = 0.

Since :math:`\lambda_t` is just :math:`\dot \lambda`, we have the boundary
condition

.. math::

   (\lambda_T^* F_{\dot{y}} )\big\vert_{t=T}  = -
     \left[
       \lambda^*(T,T)\frac{\mathrm dF_{\dot{y}}}{\mathrm dt} +
       \dot{\lambda}^* F_{\dot{y}}
     \right] \bigg\vert_{t=T}.

For the index-one DAE case, the above relation and :eq:`IDAS_adj_eqns` yield

.. math::
   (\lambda_T^* F_{\dot{y}} )\bigg\vert_{t=T} = \left[g_y - \lambda^*F_y\right]\bigg\vert_{t=T}.

For the regular implicit ODE case, :math:`F_{\dot{y}}` is invertible; thus we
have :math:`\lambda(T, T) = 0`, which leads to :math:`\lambda_T(T) = -
\dot{\lambda}(T)`. As with the final conditions for :math:`\lambda(T)` in
:eq:`IDAS_adj_eqns`, the above selection for :math:`\lambda_T(T)` is not sufficient
for index-two Hessenberg DAEs (see :cite:p:`CLPS:03` for details).

.. _IDAS.Mathematics.ASA.Checkpointing:

Checkpointing scheme
--------------------

During the backward integration, the evaluation of the right-hand side of the
adjoint system requires, at the current time, the states :math:`y` which were
computed during the forward integration phase. Since IDAS implements
variable-step integration formulas, it is unlikely that the states will be
available at the desired time and so some form of interpolation is needed. The
IDAS implementation being also variable-order, it is possible that during the
forward integration phase the order may be reduced as low as first order, which
means that there may be points in time where only :math:`y` and :math:`{\dot y}`
are available. These requirements therefore limit the choices for possible
interpolation schemes. IDAS implements two interpolation methods: a cubic
Hermite interpolation algorithm and a variable-degree polynomial interpolation
method which attempts to mimic the BDF interpolant for the forward integration.

However, especially for large-scale problems and long integration intervals, the
number and size of the vectors :math:`y` and :math:`{\dot y}` that would need to
be stored make this approach computationally intractable. Thus, IDAS settles for
a compromise between storage space and execution time by implementing a
so-called *checkpointing scheme*. At the cost of at most one additional forward
integration, this approach offers the best possible estimate of memory
requirements for adjoint sensitivity analysis. To begin with, based on the
problem size :math:`N` and the available memory, the user decides on the number
:math:`N_d` of data pairs (:math:`y`, :math:`{\dot y}`) if cubic Hermite
interpolation is selected, or on the number :math:`N_d` of :math:`y` vectors in
the case of variable-degree polynomial interpolation, that can be kept in memory
for the purpose of interpolation. Then, during the first forward integration
stage, after every :math:`N_d` integration steps a checkpoint is formed by
saving enough information (either in memory or on disk) to allow for a hot
restart, that is a restart which will exactly reproduce the forward integration.
In order to avoid storing Jacobian-related data at each checkpoint, a
reevaluation of the iteration matrix is forced before each checkpoint. At the
end of this stage, we are left with :math:`N_c` checkpoints, including one at
:math:`t_0`. During the backward integration stage, the adjoint variables are
integrated backwards from :math:`T` to :math:`t_0`, going from one checkpoint to
the previous one. The backward integration from checkpoint :math:`i+1` to
checkpoint :math:`i` is preceded by a forward integration from :math:`i` to
:math:`i+1` during which the :math:`N_d` vectors :math:`y` (and, if necessary
:math:`{\dot y}`) are generated and stored in memory for interpolation.

.. note::

   The degree of the interpolation polynomial is always that of the current BDF
   order for the forward interpolation at the first point to the right of the
   time at which the interpolated value is sought (unless too close to the
   :math:`i`-th checkpoint, in which case it uses the BDF order at the
   right-most relevant point). However, because of the FLC BDF implementation
   (see :numref:`IDAS.Mathematics.ivp_sol`), the resulting interpolation
   polynomial is only an approximation to the underlying BDF interpolant.

   The Hermite cubic interpolation option is present because it was implemented
   chronologically first and it is also used by other adjoint solvers
   (e.g. ``DASPKADJOINT``). The variable-degree polynomial is more
   memory-efficient (it requires only half of the memory storage of the cubic
   Hermite interpolation) and is more accurate.


.. figure:: /figs/idas/ckpnt.png
   :alt: Illustration of the checkpointing algorithm for generation of the
         forward solution during the integration of the adjoint system.
   :align: center

   Illustration of the checkpointing algorithm for generation of the forward
   solution during the integration of the adjoint system.

This approach transfers the uncertainty in the number of integration steps in
the forward integration phase to uncertainty in the final number of checkpoints.
However, :math:`N_c` is much smaller than the number of steps taken during the
forward integration, and there is no major penalty for writing/reading the
checkpoint data to/from a temporary file. Note that, at the end of the first
forward integration stage, interpolation data are available from the last
checkpoint to the end of the interval of integration. If no checkpoints are
necessary (:math:`N_d` is larger than the number of integration steps taken in
the solution of :eq:`IDAS_DAE_p`), the total cost of an adjoint sensitivity
computation can be as low as one forward plus one backward integration. In
addition, IDAS provides the capability of reusing a set of checkpoints for
multiple backward integrations, thus allowing for efficient computation of
gradients of several functionals :eq:`IDAS_G`.

Finally, we note that the adjoint sensitivity module in IDAS provides the
necessary infrastructure to integrate backwards in time any DAE terminal value
problem dependent on the solution of the IVP :eq:`IDAS_DAE_p`, including adjoint
systems :eq:`IDAS_adj_eqns` or :eq:`IDAS_adj1_eqns`, as well as any other quadrature ODEs
that may be needed in evaluating the integrals in :eq:`IDAS_dGdp`. In particular, for
DAE systems arising from semi-discretization of time-dependent PDEs, this
feature allows for integration of either the discretized adjoint PDE system or
the adjoint of the discretized PDE.

.. _IDAS.Mathematics.Hessian:

Second-order sensitivity analysis
=================================

In some applications (e.g., dynamically-constrained optimization) it may be
desirable to compute second-order derivative information. Considering the DAE
problem :eq:`IDAS_DAE_p` and some model output functional :math:`g(y)`, the
Hessian :math:`\mathrm d^2g/\mathrm dp^2` can be obtained in a forward
sensitivity analysis setting as

.. math:: \frac{\mathrm d^2 g}{\mathrm d p^2} = \left(g_y \otimes I_{N_p} \right ) y_{pp} + y_p^T g_{yy} y_p \, ,

where :math:`\otimes` is the Kronecker product. The second-order sensitivities
are solution of the matrix DAE system:

.. math::

   \begin{split}
     & \left( F_{\dot y} \otimes I_{N_p} \right) \cdot \dot y_{pp}  +
     \left( F_y        \otimes I_{N_p} \right) \cdot y_{pp}       +
     \left( I_N \otimes {\dot y}_p^T \right) \cdot \left( F_{\dot y \dot y} \dot y_p + F_{y \dot y} y_p \right) +
     \left( I_N \otimes y_p^T        \right) \cdot \left( F_{y \dot y}      \dot y_p + F_{y y}      y_p \right) = 0 \\
     & y_{pp}(t_0) = \frac{\partial^2 y_0}{\partial p^2} \, , \quad
     \dot y_{pp}(t_0) = \frac{\partial^2 \dot y_0}{\partial p^2} \, ,
     \end{split}

where :math:`y_p` denotes the first-order sensitivity matrix, the solution of
:math:`N_p` systems :eq:`IDAS_sens_eqns`, and :math:`y_{pp}` is a third-order tensor.
It is easy to see that, except for situations in which the number of parameters
:math:`N_p` is very small, the computational cost of this so-called
*forward-over-forward* approach is exorbitant as it requires the solution of
:math:`N_p + N_p^2` additional DAE systems of the same dimension as :eq:`IDAS_DAE_p`.

.. note::

   For the sake of simplifity in presentation, we do not include explicit
   dependencies of :math:`g` on time :math:`t` or parameters
   :math:`p`. Moreover, we only consider the case in which the dependency of the
   original DAE :eq:`IDAS_DAE_p` on the parameters :math:`p` is through its
   initial conditions only. For details on the derivation in the general case,
   see :cite:p:`OzBa:05`.

A much more efficient alternative is to compute Hessian-vector products using a
so-called *forward-over-adjoint* approach. This method is based on using the
same “trick” as the one used in computing gradients of pointwise functionals
with the adjoint method, namely applying a formal directional forward derivation
to the gradient of :eq:`IDAS_dGdp` (or the equivalent one for a pointwise functional
:math:`g(T, y(T))`). With that, the cost of computing a full Hessian is roughly
equivalent to the cost of computing the gradient with forward sensitivity
analysis. However, Hessian-vector products can be cheaply computed with one
additional adjoint solve.

As an illustration, consider the ODE problem (the derivation for the general DAE
case is too involved for the purposes of this discussion)

.. math:: {\dot y}  = f(t,\,y) \, , \quad y(t_0)  = y_0(p) \, ,

depending on some parameters :math:`p` through the initial conditions only and
consider the model functional output :math:`G(p) = \int_{t_0}^{t_f} g(t,y) \,
\mathrm dt`. It can be shown that the product between the Hessian of :math:`G` (with
respect to the parameters :math:`p`) and some vector :math:`u` can be computed
as

.. math::

   \frac{\partial^2 G}{\partial p^2} u =
     \left[ \left(\lambda^T \otimes I_{N_p} \right) y_{pp}u + y_p^T \mu \right]_{t=t_0} \, ,

where :math:`\lambda` and :math:`\mu` are solutions of

.. math::

   \begin{split}
       &-\dot\mu = f_y^T\mu + \left(\lambda^T \otimes I_n \right) f_{yy} s \, ; \quad \mu(t_f) = 0 \\
       &-\dot\lambda = f_y^T\lambda + g_y^T \, ; \quad \lambda(t_f) = 0 \\
       &\dot s = f_y s \, ; \quad s(t_0) = y_{0p} u .
     \end{split}

In the above equation, :math:`s = y_p u` is a linear combination of the columns
of the sensitivity matrix :math:`y_p`. The *forward-over-adjoint* approach
hinges crucially on the fact that :math:`s` can be computed at the cost of a
forward sensitivity analysis with respect to a single parameter (the last ODE
problem above) which is possible due to the linearity of the forward sensitivity
equations :eq:`IDAS_sens_eqns`.

Therefore (and this is also valid for the DAE case), the cost of computing the
Hessian-vector product is roughly that of two forward and two backward
integrations of a system of DAEs of size :math:`N`. For more details, including
the corresponding formulas for a pointwise model functional output, see the work
by Ozyurt and Barton :cite:p:`OzBa:05` who discuss this problem for ODE initial
value problems. As far as we know, there is no published equivalent work on DAE
problems. However, the derivations given in :cite:p:`OzBa:05` for ODE problems
can be extended to DAEs with some careful consideration given to the derivation
of proper final conditions on the adjoint systems, following the ideas presented
in :cite:p:`CLPS:03`.

To allow the *foward-over-adjoint* approach described above, IDAS provides
support for:

-  the integration of multiple backward problems depending on the same
   underlying forward problem :eq:`IDAS_DAE_p`, and

-  the integration of backward problems and computation of backward quadratures
   depending on both the states :math:`y` and forward sensitivities (for this
   particular application, :math:`s`) of the original problem :eq:`IDAS_DAE_p`.
