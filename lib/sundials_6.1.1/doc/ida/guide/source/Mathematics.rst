.. ----------------------------------------------------------------
   SUNDIALS Copyright Start
   Copyright (c) 2002-2022, Lawrence Livermore National Security
   and Southern Methodist University.
   All rights reserved.

   See the top-level LICENSE and NOTICE files for details.

   SPDX-License-Identifier: BSD-3-Clause
   SUNDIALS Copyright End
   ----------------------------------------------------------------

.. _IDA.Mathematics:

Mathematical Considerations
===========================

IDA solves the initial-value problem (IVP) for a DAE system of the general form

.. math::
   F(t,y,\dot{y}) = 0 \, ,
   \quad y(t_0) = y_0 \, , \quad \dot{y}(t_0) = \dot{y}_0 \,
   :label: IDA_DAE

where :math:`y`, :math:`\dot{y}`, and :math:`F` are vectors in :math:`{\bf
R}^N`, :math:`t` is the independent variable, :math:`\dot{y} = \mathrm dy/\mathrm dt`, and
initial values :math:`y_0`, :math:`\dot{y}_0` are given. (Often :math:`t` is
time, but it certainly need not be.)

.. _IDA.Mathematics.ivp_sol:

IVP solution
------------

Prior to integrating a DAE initial-value problem, an important requirement is
that the pair of vectors :math:`y_0` and :math:`\dot{y}_0` are both initialized
to satisfy the DAE residual :math:`F(t_0,y_0, \dot{y}_0) = 0`.  For a class of
problems that includes so-called semi-explicit index-one systems, IDA provides a
routine that computes consistent initial conditions from a user’s initial guess
:cite:p:`BHP:98`.  For this, the user must identify sub-vectors of :math:`y`
(not necessarily contiguous), denoted :math:`y_d` and :math:`y_a`, which are its
differential and algebraic parts, respectively, such that :math:`F` depends on
:math:`\dot{y}_d` but not on any components of :math:`\dot{y}_a`. The assumption
that the system is “index one” means that for a given :math:`t` and :math:`y_d`,
the system :math:`F(t,y,\dot{y}) = 0` defines :math:`y_a` uniquely. In this
case, a solver within IDA computes :math:`y_a` and :math:`\dot{y}_d` at :math:`t
= t_0`, given :math:`y_d` and an initial guess for :math:`y_a`. A second
available option with this solver also computes all of :math:`y(t_0)` given
:math:`\dot{y}(t_0)`; this is intended mainly for quasi-steady-state problems,
where :math:`\dot{y}(t_0) = 0` is given.  In both cases, IDA solves the system
:math:`F(t_0,y_0, \dot{y}_0) = 0` for the unknown components of :math:`y_0` and
:math:`\dot{y}_0`, using Newton iteration augmented with a line search global
strategy. In doing this, it makes use of the existing machinery that is to be
used for solving the linear systems during the integration, in combination with
certain tricks involving the step size (which is set artificially for this
calculation).  For problems that do not fall into either of these categories,
the user is responsible for passing consistent values, or risks failure in the
numerical integration.

The integration method used in IDA is the variable-order, variable-coefficient
BDF (Backward Differentiation Formula), in fixed-leading-coefficient form
:cite:p:`BCP:96`.  The method order ranges from 1 to 5, with the BDF of order
:math:`q` given by the multistep formula

.. math::
   \sum_{i=0}^q \alpha_{n,i}y_{n-i} = h_n \dot{y}_n \, ,
   :label: IDA_BDF

where :math:`y_n` and :math:`\dot{y}_n` are the computed approximations to
:math:`y(t_n)` and :math:`\dot{y}(t_n)`, respectively, and the step size is
:math:`h_n = t_n - t_{n-1}`.  The coefficients :math:`\alpha_{n,i}` are uniquely
determined by the order :math:`q`, and the history of the step sizes. The
application of the BDF :eq:`IDA_BDF` to the DAE system :eq:`IDA_DAE` results in a
nonlinear algebraic system to be solved at each step:

.. math::
   G(y_n) \equiv
   F \left( t_n , \, y_n , \,
      h_n^{-1} \sum_{i=0}^q \alpha_{n,i}y_{n-i} \right) = 0 \, .
   :label: IDA_DAE_nls

By default IDA solves :eq:`IDA_DAE_nls` with a Newton iteration but IDA also allows
for user-defined nonlinear solvers (see Chapter :numref:`SUNNonlinSol`). Each
Newton iteration requires the solution of a linear system of the form

.. math::
   J [y_{n(m+1)} - y_{n(m)}] = -G(y_{n(m)})  \, ,
   :label: IDA_DAE_Newtoncorr

where :math:`y_{n(m)}` is the :math:`m`-th approximation to :math:`y_n`.  Here
:math:`J` is some approximation to the system Jacobian

.. math::
   J = \frac{\partial G}{\partial y}
   = \frac{\partial F}{\partial y} +
   \alpha\frac{\partial F}{\partial \dot{y}} \, ,
   :label: IDA_DAE_Jacobian

where :math:`\alpha = \alpha_{n,0}/h_n`. The scalar :math:`\alpha` changes
whenever the step size or method order changes.

For the solution of the linear systems within the Newton iteration, IDA provides
several choices, including the option of a user-supplied linear solver (see
Chapter :numref:`SUNLinSol`). The linear solvers distributed with SUNDIALS are
organized in two families, a *direct* family comprising direct linear solvers
for dense, banded, or sparse matrices and a *spils* family comprising scaled
preconditioned iterative (Krylov) linear solvers.  The methods offered through
these modules are as follows:

* dense direct solvers, using either an internal implementation or
  a BLAS/LAPACK implementation (serial or threaded vector modules only),

* band direct solvers, using either an internal implementation or
  a BLAS/LAPACK implementation (serial or threaded vector modules only),

* sparse direct solver interfaces, using either the KLU sparse solver
  library :cite:p:`DaPa:10,KLU_site`, or the thread*enabled SuperLU_MT sparse
  solver library :cite:p:`Li:05,DGL:99,SuperLUMT_site` (serial or threaded
  vector modules only) [Note that users will need to download and install the
  KLU or SuperLU_MT packages independent of IDA],

* SPGMR, a scaled preconditioned GMRES (Generalized Minimal Residual method)
  solver with or without restarts,

* SPFGMR, a scaled preconditioned FGMRES (Flexible Generalized
  Minimal Residual method) solver with or without restarts,

* SPBCG, a scaled preconditioned Bi-CGStab (Bi-Conjugate Gradient Stable
  method) solver,

* SPTFQMR, a scaled preconditioned TFQMR (Transpose-Free Quasi-Minimal
  Residual method) solver, or

* PCG, a scaled preconditioned CG (Conjugate Gradient method) solver.

For large stiff systems, where direct methods are not feasible, the combination
of a BDF integrator and a preconditioned Krylov method yields a powerful tool
because it combines established methods for stiff integration, nonlinear
iteration, and Krylov (linear) iteration with a problem-specific treatment of
the dominant source of stiffness, in the form of the user-supplied
preconditioner matrix :cite:p:`BrHi:89`.  For the *spils* linear solvers with
IDA, preconditioning is allowed only on the left (see
:numref:`IDA.Mathematics.preconditioning`).  Note that the dense, band, and sparse
direct linear solvers can only be used with serial and threaded vector
representations.

In the process of controlling errors at various levels, IDA uses a weighted
root-mean-square norm, denoted :math:`\|\cdot\|_{\mbox{WRMS}}`, for all
error-like quantities. The multiplicative weights used are based on the current
solution and on the relative and absolute tolerances input by the user, namely

.. math::
   W_i = \frac{1}{\text{rtol} \cdot |y_i| + \text{atol}_i }\, .
   :label: IDA_errwt

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

The default stopping test for nonlinear solver iterations in IDA ensures that
the iteration error :math:`y_n - y_{n(m)}` is small relative to :math:`y`
itself. For this, we estimate the linear convergence rate at all iterations
:math:`m>1` as

.. math:: R = \left( \frac{\delta_m}{\delta_1} \right)^{\frac{1}{m-1}} \, ,

where the :math:`\delta_m = y_{n(m)} - y_{n(m-1)}` is the correction at
iteration :math:`m=1,2,\ldots`. The nonlinear solver iteration is halted if
:math:`R>0.9`.  The convergence test at the :math:`m`-th iteration is then

.. math::
   S \| \delta_m \| < 0.33 \, ,
   :label: IDA_DAE_nls_test

where :math:`S = R/(R-1)` whenever :math:`m>1` and :math:`R\le 0.9`. The user
has the option of changing the constant in the convergence test from its default
value of :math:`0.33`.  The quantity :math:`S` is set to :math:`S=20` initially
and whenever :math:`J` or :math:`P` is updated, and it is reset to :math:`S=100`
on a step with :math:`\alpha \neq \bar\alpha`.  Note that at :math:`m=1`, the
convergence test :eq:`IDA_DAE_nls_test` uses an old value for :math:`S`. Therefore,
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
or :ref:`SUNMATRIX_BAND <SUNMatrix.Band>` matrix objects,
the Jacobian :math:`J` defined in :eq:`IDA_DAE_Jacobian` can be either supplied by
the user or IDA can compute :math:`J` internally by difference quotients. In the
latter case, we use the approximation

.. math::
   \begin{gathered}
     J_{ij} = [F_i(t,y+\sigma_j e_j,\dot{y}+\alpha\sigma_j e_j) -
               F_i(t,y,\dot{y})]/\sigma_j \, , \text{ with}\\
     \sigma_j = \sqrt{U} \max \left\{ |y_j|, |h\dot{y}_j|,1/W_j \right\}
                \mbox{sign}(h \dot{y}_j) \, ,\end{gathered}

where :math:`U` is the unit roundoff, :math:`h` is the current step size, and
:math:`W_j` is the error weight for the component :math:`y_j` defined by
:eq:`IDA_errwt`.  We note that with sparse and user-supplied matrix objects,
the Jacobian *must* be supplied by a user routine.

In the case of an iterative linear solver, if a routine for :math:`Jv` is not
supplied, such products are approximated by

.. math:: Jv = [F(t,y+\sigma v,\dot{y}+\alpha\sigma v) - F(t,y,\dot{y})]/\sigma \, ,

where the increment :math:`\sigma = 1/\|v\|`. As an option, the user can specify
a constant factor that is inserted into this expression for :math:`\sigma`.

During the course of integrating the system, IDA computes an estimate of the
local truncation error, LTE, at the :math:`n`-th time step, and requires this to
satisfy the inequality

.. math:: \| \mbox{LTE} \|_{\mbox{WRMS}} \leq 1 \, .

Asymptotically, LTE varies as :math:`h^{q+1}` at step size :math:`h` and order
:math:`q`, as does the predictor-corrector difference :math:`\Delta_n \equiv
y_n-y_{n(0)}`.  Thus there is a constant :math:`C` such that

.. math:: \mbox{LTE} = C \Delta_n + O(h^{q+2}) \, ,

and so the norm of LTE is estimated as :math:`|C| \cdot \|\Delta_n\|`.  In
addition, IDA requires that the error in the associated polynomial interpolant
over the current step be bounded by 1 in norm. The leading term of the norm of
this error is bounded by :math:`\bar{C} \|\Delta_n\|` for another constant
:math:`\bar{C}`. Thus the local error test in IDA is

.. math::
   \max\{ |C|, \bar{C} \} \|\Delta_n\| \leq 1 \, .
   :label: IDA_lerrtest

A user option is available by which the algebraic components of the error vector
are omitted from the test :eq:`IDA_lerrtest`, if these have been so identified.

In IDA, the local error test is tightly coupled with the logic for selecting the
step size and order. First, there is an initial phase that is treated specially;
for the first few steps, the step size is doubled and the order raised (from its
initial value of 1) on every step, until (a) the local error test :eq:`IDA_lerrtest`
fails, (b) the order is reduced (by the rules given below), or (c) the order
reaches 5 (the maximum). For step and order selection on the general step, IDA
uses a different set of local error estimates, based on the asymptotic behavior
of the local error in the case of fixed step sizes.  At each of the orders
:math:`q'` equal to :math:`q`, :math:`q-1` (if :math:`q > 1`), :math:`q-2` (if
:math:`q > 2`), or :math:`q+1` (if :math:`q < 5`), there are constants
:math:`C(q')` such that the norm of the local truncation error at order
:math:`q'` satisfies

.. math:: \mbox{LTE}(q') = C(q') \| \phi(q'+1) \| + O(h^{q'+2}) \, ,

where :math:`\phi(k)` is a modified divided difference of order :math:`k` that
is retained by IDA (and behaves asymptotically as :math:`h^k`).  Thus the local
truncation errors are estimated as ELTE\ :math:`(q') = C(q')\|\phi(q'+1)\|` to
select step sizes. But the choice of order in IDA is based on the requirement
that the scaled derivative norms, :math:`\|h^k y^{(k)}\|`, are monotonically
decreasing with :math:`k`, for :math:`k` near :math:`q`. These norms are again
estimated using the :math:`\phi(k)`, and in fact

.. math:: \|h^{q'+1} y^{(q'+1)}\| \approx T(q') \equiv (q'+1) \mbox{ELTE}(q') \, .

The step/order selection begins with a test for monotonicity that is made even
*before* the local error test is performed. Namely, the order is reset to
:math:`q' = q-1` if (a) :math:`q=2` and :math:`T(1)\leq T(2)/2`, or (b) :math:`q
> 2` and :math:`\max\{T(q-1),T(q-2)\} \leq T(q)`; otherwise :math:`q' = q`. Next
the local error test :eq:`IDA_lerrtest` is performed, and if it fails, the step is
redone at order :math:`q\leftarrow q'` and a new step size :math:`h'`. The
latter is based on the :math:`h^{q+1}` asymptotic behavior of
:math:`\mbox{ELTE}(q)`, and, with safety factors, is given by

.. math:: \eta = h'/h = 0.9/[2 \, \mbox{ELTE}(q)]^{1/(q+1)} \, .

The value of :math:`\eta` is adjusted so that :math:`0.25 \leq \eta \leq 0.9`
before setting :math:`h \leftarrow h' = \eta h`. If the local error test fails a
second time, IDA uses :math:`\eta = 0.25`, and on the third and subsequent
failures it uses :math:`q = 1` and :math:`\eta = 0.25`. After 10 failures, IDA
returns with a give-up message.

As soon as the local error test has passed, the step and order for the next step
may be adjusted. No such change is made if :math:`q' = q-1` from the prior test,
if :math:`q = 5`, or if :math:`q` was increased on the previous step. Otherwise,
if the last :math:`q+1` steps were taken at a constant order :math:`q < 5` and a
constant step size, IDA considers raising the order to :math:`q+1`. The logic is
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

IDA permits the user to impose optional inequality constraints on individual
components of the solution vector :math:`y`. Any of the following four
constraints can be imposed: :math:`y_i > 0`, :math:`y_i < 0`, :math:`y_i \geq
0`, or :math:`y_i \leq 0`.  The constraint satisfaction is tested after a
successful nonlinear system solution.  If any constraint fails, we declare a
convergence failure of the nonlinear iteration and reduce the step size. Rather
than cutting the step size by some arbitrary factor, IDA estimates a new step
size :math:`h'` using a linear approximation of the components in :math:`y` that
failed the constraint test (including a safety factor of :math:`0.9` to cover
the strict inequality case). These additional constraints are also imposed
during the calculation of consistent initial conditions.  If a step fails to
satisfy the constraints repeatedly within a step attempt then the integration is
halted and an error is returned. In this case the user may need to employ other
strategies as discussed in :numref:`IDA.Usage.CC.callable_fct_sim.idatolerances` to
satisfy the inequality constraints.

Normally, IDA takes steps until a user-defined output value :math:`t =
t_{\mbox{out}}` is overtaken, and then computes :math:`y(t_{\mbox{out}})` by
interpolation. However, a “one step” mode option is available, where control
returns to the calling program after each step. There are also options to force
IDA not to integrate past a given stopping point :math:`t = t_{\mbox{stop}}`.

.. _IDA.Mathematics.preconditioning:

Preconditioning
---------------

When using a nonlinear solver that requires the solution of a linear system of
the form :math:`J \Delta y = - G` (e.g., the default Newton iteration), IDA
makes repeated use of a linear solver.  If this linear system solve is done with
one of the scaled preconditioned iterative linear solvers supplied with
SUNDIALS, these solvers are rarely successful if used without preconditioning;
it is generally necessary to precondition the system in order to obtain
acceptable efficiency.  A system :math:`A x = b` can be preconditioned on the
left, on the right, or on both sides. The Krylov method is then applied to a
system with the matrix :math:`P^{-1}A`, or :math:`AP^{-1}`, or :math:`P_L^{-1} A
P_R^{-1}`, instead of :math:`A`.  However, within IDA, preconditioning is
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

Typical preconditioners used with IDA are based on approximations to the
iteration matrix of the systems involved; in other words, :math:`P \approx
\dfrac{\partial F}{\partial y} + \alpha\dfrac{\partial F}{\partial \dot{y}}`,
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

.. _IDA.Mathematics.rootfinding:

Rootfinding
-----------

The IDA solver has been augmented to include a rootfinding feature. This means
that, while integratnuming the Initial Value Problem :eq:`IDA_DAE`, IDA can also
find the roots of a set of user-defined functions :math:`g_i(t,y,\dot{y})` that
depend on :math:`t`, the solution vector :math:`y = y(t)`, and its :math:`t-`\
derivative :math:`\dot{y}(t)`. The number of these root functions is arbitrary,
and if more than one :math:`g_i` is found to have a root in any given interval,
the various root locations are found and reported in the order that they occur
on the :math:`t` axis, in the direction of integration.

Generally, this rootfinding feature finds only roots of odd multiplicity,
corresponding to changes in sign of :math:`g_i(t,y(t),\dot{y}(t))`, denoted
:math:`g_i(t)` for short. If a user root function has a root of even
multiplicity (no sign change), it will probably be missed by IDA. If such a root
is desired, the user should reformulate the root function so that it changes
sign at the desired root.

The basic scheme used is to check for sign changes of any :math:`g_i(t)` over
each time step taken, and then (when a sign change is found) to home in on the
root (or roots) with a modified secant method :cite:p:`HeSh:80`.  In addition,
each time :math:`g` is computed, IDA checks to see if :math:`g_i(t) = 0`
exactly, and if so it reports this as a root. However, if an exact zero of any
:math:`g_i` is found at a point :math:`t`, IDA computes :math:`g` at :math:`t +
\delta` for a small increment :math:`\delta`, slightly further in the direction
of integration, and if any :math:`g_i(t + \delta)=0` also, IDA stops and reports
an error. This way, each time IDA takes a time step, it is guaranteed that the
values of all :math:`g_i` are nonzero at some past value of :math:`t`, beyond
which a search for roots is to be done.

At any given time in the course of the time-stepping, after suitable checking
and adjusting has been done, IDA has an interval :math:`(t_{lo},t_{hi}]` in
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
