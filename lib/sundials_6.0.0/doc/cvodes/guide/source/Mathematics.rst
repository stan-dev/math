.. ----------------------------------------------------------------
   SUNDIALS Copyright Start
   Copyright (c) 2002-2021, Lawrence Livermore National Security
   and Southern Methodist University.
   All rights reserved.

   See the top-level LICENSE and NOTICE files for details.

   SPDX-License-Identifier: BSD-3-Clause
   SUNDIALS Copyright End
   ----------------------------------------------------------------

.. _CVODES.Mathematics:

***************************
Mathematical Considerations
***************************

CVODES solves ODE initial value problems (IVPs) in real
:math:`N`-space, which we write in the abstract form

.. math::
   \dot{y} = f(t,y) \, ,\quad y(t_0) = y_0 \,
   :label: CVODES_ivp

where :math:`y \in \mathbb{R}^N` and
:math:`f: \mathbb{R} \times \mathbb{R}^N \rightarrow \mathbb{R}^N`.
Here we use :math:`\dot{y}` to denote :math:`\mathrm dy/\mathrm dt`. While we use
:math:`t` to denote the independent variable, and usually this is time,
it certainly need not be. CVODES solves both stiff and nonstiff
systems. Roughly speaking, stiffness is characterized by the presence of
at least one rapidly damped mode, whose time constant is small compared
to the time scale of the solution itself.

Additionally, if :eq:`CVODES_ivp` depends on some parameters :math:`p \in {\bf R}^{N_p}`, i.e.

.. math::
   \begin{split}
   &{\dot{y}}  = f(t,\,y,\,p) \\
   &y(t_0)  = y_0(p) \, ,
   \end{split}
   :label: CVODES_ivp_p

CVODES can also compute first order derivative information, performing either
*forward sensitivity analysis* or *adjoint sensitivity analysis*. In the first
case, CVODES computes the sensitivities of the solution with respect to the
parameters :math:`p`, while in the second case, CVODES computes the gradient of
a *derived function* with respect to the parameters :math:`p`.


.. _CVODES.Mathematics.ivp_sol:

IVP solution
============

The methods used in CVODES are variable-order, variable-step
multistep methods, based on formulas of the form

.. math::
   \sum_{i = 0}^{K_1} \alpha_{n,i} y^{n-i} +
        h_n \sum_{i = 0}^{K_2} \beta_{n,i} {\dot{y}}^{n-i} = 0 \, .
   :label: CVODES_lmm

Here the :math:`y^n` are computed approximations to :math:`y(t_n)`, and
:math:`h_n = t_n - t_{n-1}` is the step size. The user of CVODES must
choose appropriately one of two multistep methods. For nonstiff
problems, CVODES includes the Adams-Moulton formulas, characterized
by :math:`K_1 = 1` and :math:`K_2 = q-1` above, where the order
:math:`q` varies between :math:`1` and :math:`12`. For stiff problems,
CVODES includes the Backward Differentiation Formulas (BDF) in
so-called fixed-leading coefficient (FLC) form, given by :math:`K_1 = q`
and :math:`K_2 = 0`, with order :math:`q` varying between :math:`1` and
:math:`5`. The coefficients are uniquely determined by the method type,
its order, the recent history of the step sizes, and the normalization
:math:`\alpha_{n,0} = -1`. See :cite:p:`ByHi:75` and
:cite:p:`JaSD:80`.

For either choice of formula, a nonlinear system must be solved
(approximately) at each integration step. This nonlinear system can be
formulated as either a rootfinding problem

.. math::
   F(y^n) \equiv y^n - h_n \beta_{n,0} f(t_n,y^n) - a_n = 0 \, ,
   :label: CVODES_nonlinear

or as a fixed-point problem

.. math::
   G(y^n) \equiv h_n \beta_{n,0} f(t_n,y^n) + a_n = y^n \, .
   :label: CVODES_nonlinear_fixedpoint

where
:math:`a_n\equiv\sum_{i>0}(\alpha_{n,i}y^{n-i}+h_n\beta_{n,i} {\dot{y}}^{n-i})`.
CVODES provides several nonlinear solver choices as well as the
option of using a user-defined nonlinear solver (see
:numref:`SUNNonlinSol`). By default CVODES solves :eq:`CVODES_nonlinear` with a
*Newton iteration* which requires the solution of linear systems

.. math::
   M [y^{n(m+1)} - y^{n(m)}] = -F(y^{n(m)}) \, ,
   :label: CVODES_Newton

in which

.. math::
   M \approx I - \gamma J \, ,
   \quad J = \partial f / \partial y \, ,
   \quad \mbox{and} \quad
   \gamma = h_n \beta_{n,0} \, .
   :label: CVODES_Newtonmat

The exact variation of the Newton iteration depends on the choice of linear
solver and is discussed below and in :numref:`SUNNonlinSol.Newton`. For nonstiff
systems, a *fixed-point iteration* (previously referred to as a functional
iteration in this guide) solving :eq:`CVODES_nonlinear_fixedpoint` is also
available. This involves evaluations of :math:`f` only and can (optionally) use
Anderson’s method :cite:p:`Anderson65, Walker-Ni09, Fang-Saad09, LWWY11` to
accelerate convergence (see :numref:`SUNNonlinSol.FixedPoint` for more details).
For any nonlinear solver, the initial guess for the iteration is a predicted
value :math:`y^{n(0)}` computed explicitly from the available history data.

For nonlinear solvers that require the solution of the linear system
:eq:`CVODES_Newton` (e.g., the default Newton iteration),
CVODES provides several linear solver choices, including the option
of a user-supplied linear solver module (see
:numref:`SUNLinSol`). The linear solver modules distributed
with SUNDIALS are organized in two families, a *direct* family
comprising direct linear solvers for dense, banded, or sparse matrices,
and a *spils* family comprising scaled preconditioned iterative (Krylov)
linear solvers. The methods offered through these modules are as
follows:

* dense direct solvers, including an internal implementation, an interface to
  BLAS/LAPACK, an interface to MAGMA :cite:p:`magma_ref` and an interface to
  the oneMKL library :cite:p:`oneAPI_site`,

* band direct solvers, including an internal implementation or an interface to BLAS/LAPACK,

* sparse direct solver interfaces to various libraries, including KLU :cite:p:`DaPa:10, KLU_site`,
  SuperLU_MT :cite:p:`Li:05,DGL:99,SuperLUMT_site`, SuperLU_Dist
  :cite:p:`GDL:07,LD:03,SLUUG:99,SuperLUDIST_site`, and cuSPARSE :cite:p:`cuSPARSE_site`,

* SPGMR, a scaled preconditioned GMRES (Generalized Minimal Residual method) solver,

* SPFGMR, a scaled preconditioned FGMRES (Flexible Generalized Minimal Residual method) solver,

* SPBCG, a scaled preconditioned Bi-CGStab (Bi-Conjugate Gradient Stable method) solver,

* SPTFQMR, a scaled preconditioned TFQMR (Transpose-Free Quasi-Minimal Residual method) solver, or

* PCG, a scaled preconditioned CG (Conjugate Gradient method) solver.

For large stiff systems, where direct methods are often not feasible,
the combination of a BDF integrator and a preconditioned Krylov method
yields a powerful tool because it combines established methods for stiff
integration, nonlinear iteration, and Krylov (linear) iteration with a
problem-specific treatment of the dominant source of stiffness, in the
form of the user-supplied preconditioner matrix
:cite:p:`BrHi:89`.

In addition, CVODES also provides a linear solver module which only
uses a diagonal approximation of the Jacobian matrix.

In the process of controlling errors at various levels, CVODES uses a
weighted root-mean-square norm, denoted
:math:`|\cdot|_{\text{WRMS}}`, for all error-like
quantities. The multiplicative weights used are based on the current
solution and on the relative and absolute tolerances input by the user,
namely

.. math::
   W_i = 1 / [\text{rtol} \cdot |y_i| + \text{atol}_i ] \, .
   :label: CVODES_errwt

Because :math:`1/W_i` represents a tolerance in the component
:math:`y_i`, a vector whose norm is 1 is regarded as “small.” For
brevity, we will usually drop the subscript WRMS on norms in what
follows.

In the case of a matrix-based linear solver, the default Newton
iteration is a Modified Newton iteration, in that the iteration matrix
:math:`M` is fixed throughout the nonlinear iterations. However, in the
case that a matrix-free iterative linear solver is used, the default
Newton iteration is an Inexact Newton iteration, in which :math:`M` is
applied in a matrix-free manner, with matrix-vector products :math:`Jv`
obtained by either difference quotients or a user-supplied routine. With
the default Newton iteration, the matrix :math:`M` and preconditioner
matrix :math:`P` are updated as infrequently as possible to balance the
high costs of matrix operations against other costs. Specifically, this
matrix update occurs when:

   * starting the problem,
   * more than 20 steps have been taken since the last update,
   * the value :math:`\bar{\gamma}` of :math:`\gamma` at the last update satisfies :math:`|\gamma/\bar{\gamma} - 1| > 0.3`,
   * a non-fatal convergence failure just occurred, or
   * an error test failure just occurred.

When forced by a convergence failure, an update of :math:`M` or
:math:`P` may or may not involve a reevaluation of :math:`J` (in
:math:`M`) or of Jacobian data (in :math:`P`), depending on whether
Jacobian error was the likely cause of the failure. More generally, the
decision is made to reevaluate :math:`J` (or instruct the user to
reevaluate Jacobian data in :math:`P`) when:

   * starting the problem,
   * more than 50 steps have been taken since the last evaluation,
   * a convergence failure occurred with an outdated matrix, and the value :math:`\bar{\gamma}` of :math:`\gamma` at the last update satisfies :math:`|\gamma/\bar{\gamma} - 1| < 0.2`, or
   * a convergence failure occurred that forced a step size reduction.

The default stopping test for nonlinear solver iterations is related to
the subsequent local error test, with the goal of keeping the nonlinear
iteration errors from interfering with local error control. As described
below, the final computed value :math:`y^{n(m)}` will have to satisfy a
local error test :math:`\|y^{n(m)} - y^{n(0)}\| \leq \epsilon`. Letting
:math:`y^n` denote the exact solution of :eq:`CVODES_nonlinear`, we want to ensure that the iteration
error :math:`y^n - y^{n(m)}` is small relative to :math:`\epsilon`,
specifically that it is less than :math:`0.1 \epsilon`. (The safety
factor :math:`0.1` can be changed by the user.) For this, we also
estimate the linear convergence rate constant :math:`R` as follows. We
initialize :math:`R` to 1, and reset :math:`R = 1` when :math:`M` or
:math:`P` is updated. After computing a correction
:math:`\delta_m = y^{n(m)}-y^{n(m-1)}`, we update :math:`R` if
:math:`m > 1` as

.. math:: R \leftarrow \max\{0.3R , \|\delta_m\| / \|\delta_{m-1}\| \} \, .

Now we use the estimate

.. math::

   \| y^n - y^{n(m)} \| \approx \| y^{n(m+1)} - y^{n(m)} \|
     \approx R \| y^{n(m)} - y^{n(m-1)} \|  =  R \|\delta_m \| \, .

Therefore the convergence (stopping) test is

.. math:: R \|\delta_m\| < 0.1 \epsilon \, .

We allow at most 3 iterations (but this limit can be changed by the
user). We also declare the iteration diverged if any
:math:`\|\delta_m\| / \|\delta_{m-1}\| > 2` with :math:`m > 1`. If convergence fails with
:math:`J` or :math:`P` current, we are forced to reduce the step size,
and we replace :math:`h_n` by :math:`h_n/4`. The integration is halted
after a preset number of convergence failures; the default value of this
limit is 10, but this can be changed by the user.

When an iterative method is used to solve the linear system, its errors
must also be controlled, and this also involves the local error test
constant. The linear iteration error in the solution vector
:math:`\delta_m` is approximated by the preconditioned residual vector.
Thus to ensure (or attempt to ensure) that the linear iteration errors
do not interfere with the nonlinear error and local integration error
controls, we require that the norm of the preconditioned residual be
less than :math:`0.05 \cdot (0.1 \epsilon)`.

When the Jacobian is stored using either the :ref:`SUNMATRIX_DENSE <SUNMatrix.Dense>`
or :ref:`SUNMATRIX_BAND <SUNMatrix.Band>` matrix
objects, the Jacobian may be supplied by a user routine, or approximated
by difference quotients, at the user’s option. In the latter case, we
use the usual approximation

.. math:: J_{ij} = [f_i(t,y+\sigma_j e_j) - f_i(t,y)]/\sigma_j \, .

The increments :math:`\sigma_j` are given by

.. math:: \sigma_j = \max\left\{\sqrt{U} \; |y_j| , \sigma_0 / W_j \right\} \, ,

where :math:`U` is the unit roundoff, :math:`\sigma_0` is a
dimensionless value, and :math:`W_j` is the error weight defined in
:eq:`CVODES_errwt`. In the dense case, this scheme requires
:math:`N` evaluations of :math:`f`, one for each column of :math:`J`. In
the band case, the columns of :math:`J` are computed in groups, by the
Curtis-Powell-Reid algorithm, with the number of :math:`f` evaluations
equal to the bandwidth.

We note that with sparse and user-supplied ``SUNMatrix`` objects, the
Jacobian *must* be supplied by a user routine.

In the case of a Krylov method, preconditioning may be used on the left,
on the right, or both, with user-supplied routines for the
preconditioning setup and solve operations, and optionally also for the
required matrix-vector products :math:`Jv`. If a routine for :math:`Jv`
is not supplied, these products are computed as

.. math::
   Jv = [f(t,y+\sigma v) - f(t,y)]/\sigma \, .
   :label: CVODES_jacobv

The increment :math:`\sigma` is :math:`1/\|v\|`, so that
:math:`\sigma v` has norm 1.

A critical part of CVODES — making it an ODE “solver” rather than
just an ODE method, is its control of local error. At every step, the
local error is estimated and required to satisfy tolerance conditions,
and the step is redone with reduced step size whenever that error test
fails. As with any linear multistep method, the local truncation error
LTE, at order :math:`q` and step size :math:`h`, satisfies an asymptotic
relation

.. math:: \mbox{LTE} = C h^{q+1} y^{(q+1)} + O(h^{q+2})

for some constant :math:`C`, under mild assumptions on the step sizes. A
similar relation holds for the error in the predictor :math:`y^{n(0)}`.
These are combined to get a relation

.. math:: \mbox{LTE} = C' [y^n - y^{n(0)}] + O(h^{q+2}) \, .

The local error test is simply :math:`|\mbox{LTE}| \leq 1`. Using the
above, it is performed on the predictor-corrector difference
:math:`\Delta_n \equiv y^{n(m)} - y^{n(0)}` (with :math:`y^{n(m)}` the
final iterate computed), and takes the form

.. math:: \|\Delta_n\| \leq \epsilon \equiv 1/|C'| \, .

If this test passes, the step is considered successful. If it fails, the
step is rejected and a new step size :math:`h'` is computed based on the
asymptotic behavior of the local error, namely by the equation

.. math:: (h'/h)^{q+1} \|\Delta_n\| = \epsilon/6 \, .

Here 1/6 is a safety factor. A new attempt at the step is made, and the
error test repeated. If it fails three times, the order :math:`q` is
reset to 1 (if :math:`q > 1`), or the step is restarted from scratch (if
:math:`q = 1`). The ratio :math:`h'/h` is limited above to 0.2 after two
error test failures, and limited below to 0.1 after three. After seven
failures, CVODES returns to the user with a give-up message.

In addition to adjusting the step size to meet the local error test,
CVODES periodically adjusts the order, with the goal of maximizing
the step size. The integration starts out at order 1 and varies the
order dynamically after that. The basic idea is to pick the order
:math:`q` for which a polynomial of order :math:`q` best fits the
discrete data involved in the multistep method. However, if either a
convergence failure or an error test failure occurred on the step just
completed, no change in step size or order is done. At the current order
:math:`q`, selecting a new step size is done exactly as when the error
test fails, giving a tentative step size ratio

.. math:: h'/h = (\epsilon / 6 \|\Delta_n\| )^{1/(q+1)} \equiv \eta_q \, .

We consider changing order only after taking :math:`q+1` steps at order
:math:`q`, and then we consider only orders :math:`q' = q - 1` (if
:math:`q > 1`) or :math:`q' = q + 1` (if :math:`q < 5`). The local
truncation error at order :math:`q'` is estimated using the history
data. Then a tentative step size ratio is computed on the basis that
this error, LTE\ :math:`(q')`, behaves asymptotically as
:math:`h^{q'+1}`. With safety factors of 1/6 and 1/10 respectively,
these ratios are:

.. math:: h'/h = [1 / 6 \|\mbox{LTE}(q-1)\| ]^{1/q} \equiv \eta_{q-1}

and

.. math:: h'/h = [1 / 10 \|\mbox{LTE}(q+1)\| ]^{1/(q+2)} \equiv \eta_{q+1} \, .

The new order and step size are then set according to

.. math:: \eta = \max\{\eta_{q-1},\eta_q,\eta_{q+1}\} ~,~~ h' = \eta h \, ,

with :math:`q'` set to the index achieving the above maximum. However,
if we find that :math:`\eta < 1.5`, we do not bother with the change.
Also, :math:`h'/h` is always limited to 10, except on the first step,
when it is limited to :math:`10^4`.

The various algorithmic features of CVODES described above, as
inherited from VODE and VODPK, are documented in
:cite:p:`BBH:89,Byr:92,Hin:00`. They are also summarized in
:cite:p:`HBGLSSW:05`.

CVODES permits the user to impose optional inequality constraints on
individual components of the solution vector :math:`y`. Any of the
following four constraints can be imposed: :math:`y_i > 0`,
:math:`y_i < 0`, :math:`y_i \geq 0`, or :math:`y_i \leq 0`. The
constraint satisfaction is tested after a successful nonlinear system
solution. If any constraint fails, we declare a convergence failure of
the Newton iteration and reduce the step size. Rather than cutting the
step size by some arbitrary factor, CVODES estimates a new step size
:math:`h'` using a linear approximation of the components in :math:`y`
that failed the constraint test (including a safety factor of
:math:`0.9` to cover the strict inequality case). If a step fails to
satisfy the constraints repeatedly within a step attempt or fails with
the minimum step size then the integration is halted and an error is
returned. In this case the user may need to employ other strategies as
discussed in :numref:`CVODES.Usage.SIM.user_callable.cvtolerances` to satisfy
the inequality constraints.

Normally, CVODES takes steps until a user-defined output value
:math:`t = t_{\text{out}}` is overtaken, and then it
computes :math:`y(t_{\text{out}})` by interpolation.
However, a “one step” mode option is available, where control returns to
the calling program after each step. There are also options to force
CVODES not to integrate past a given stopping point
:math:`t = t_{\text{stop}}`.


.. _CVODES.Mathematics.preconditioning:

Preconditioning
===============

When using a nonlinear solver that requires the solution of the linear
system (:numref:`SUNNonlinSol.Newton`) (e.g., the default Newton
iteration), CVODES makes repeated use of a linear solver to solve
linear systems of the form :math:`M x = - r`, where :math:`x` is a
correction vector and :math:`r` is a residual vector. If this linear
system solve is done with one of the scaled preconditioned iterative
linear solvers supplied with SUNDIALS, these solvers are rarely
successful if used without preconditioning; it is generally necessary to
precondition the system in order to obtain acceptable efficiency. A
system :math:`M x = b` can be preconditioned on the left, as
:math:`(P^{-1}M) x = P^{-1} b`; on the right, as
:math:`(M P^{-1}) P x = b`; or on both sides, as
:math:`(P_L^{-1} M P_R^{-1}) P_R x = P_L^{-1}b`. The Krylov method is
then applied to a system with the matrix :math:`P^{-1}M`, or
:math:`MP^{-1}`, or :math:`P_L^{-1} M P_R^{-1}`, instead of :math:`M`.
In order to improve the convergence of the Krylov iteration, the
preconditioner matrix :math:`P`, or the product :math:`P_L P_R` in the
last case, should in some sense approximate the system matrix :math:`M`.
Yet at the same time, in order to be cost-effective, the matrix
:math:`P`, or matrices :math:`P_L` and :math:`P_R`, should be reasonably
efficient to evaluate and solve. Finding a good point in this tradeoff
between rapid convergence and low cost can be very difficult. Good
choices are often problem-dependent (for example, see
:cite:p:`BrHi:89` for an extensive study of preconditioners
for reaction-transport systems).

Most of the iterative linear solvers supplied with SUNDIALS allow
for preconditioning either side, or on both sides, although we know of
no situation where preconditioning on both sides is clearly superior to
preconditioning on one side only (with the product :math:`P_L P_R`).
Moreover, for a given preconditioner matrix, the merits of left
vs. right preconditioning are unclear in general, and the user should
experiment with both choices. Performance will differ because the
inverse of the left preconditioner is included in the linear system
residual whose norm is being tested in the Krylov algorithm. As a rule,
however, if the preconditioner is the product of two matrices, we
recommend that preconditioning be done either on the left only or the
right only, rather than using one factor on each side.

Typical preconditioners used with CVODES are based on approximations
to the system Jacobian, :math:`J = \partial f / \partial y`. Since the
matrix involved is :math:`M = I - \gamma J`, any approximation
:math:`\bar{J}` to :math:`J` yields a matrix that is of potential use as
a preconditioner, namely :math:`P = I - \gamma \bar{J}`. Because the
Krylov iteration occurs within a nonlinear solver iteration and further
also within a time integration, and since each of these iterations has
its own test for convergence, the preconditioner may use a very crude
approximation, as long as it captures the dominant numerical feature(s)
of the system. We have found that the combination of a preconditioner
with the Newton-Krylov iteration, using even a fairly poor approximation
to the Jacobian, can be surprisingly superior to using the same matrix
without Krylov acceleration (i.e., a modified Newton iteration), as well
as to using the Newton-Krylov method with no preconditioning.

.. _CVODES.Mathematics.stablimit:

BDF stability limit detection
=============================

CVODES includes an algorithm, STALD (STAbility Limit Detection),
which provides protection against potentially unstable behavior of the
BDF multistep integration methods in certain situations, as described
below.

When the BDF option is selected, CVODES uses Backward
Differentiation Formula methods of orders 1 to 5. At order 1 or 2, the
BDF method is A-stable, meaning that for any complex constant
:math:`\lambda` in the open left half-plane, the method is
unconditionally stable (for any step size) for the standard scalar model
problem :math:`\dot{y} = \lambda y`. For an ODE system, this means that,
roughly speaking, as long as all modes in the system are stable, the
method is also stable for any choice of step size, at least in the sense
of a local linear stability analysis.

At orders 3 to 5, the BDF methods are not A-stable, although they are
*stiffly stable*. In each case, in order for the method to be stable at
step size :math:`h` on the scalar model problem, the product
:math:`h\lambda` must lie within a *region of absolute stability*. That
region excludes a portion of the left half-plane that is concentrated
near the imaginary axis. The size of that region of instability grows as
the order increases from 3 to 5. What this means is that, when running
BDF at any of these orders, if an eigenvalue :math:`\lambda` of the
system lies close enough to the imaginary axis, the step sizes :math:`h`
for which the method is stable are limited (at least according to the
linear stability theory) to a set that prevents :math:`h\lambda` from
leaving the stability region. The meaning of *close enough* depends on
the order. At order 3, the unstable region is much narrower than at
order 5, so the potential for unstable behavior grows with order.

System eigenvalues that are likely to run into this instability are ones
that correspond to weakly damped oscillations. A pure undamped
oscillation corresponds to an eigenvalue on the imaginary axis. Problems
with modes of that kind call for different considerations, since the
oscillation generally must be followed by the solver, and this requires
step sizes (:math:`h \sim 1/\nu`, where :math:`\nu` is the frequency)
that are stable for BDF anyway. But for a weakly damped oscillatory
mode, the oscillation in the solution is eventually damped to the noise
level, and at that time it is important that the solver not be
restricted to step sizes on the order of :math:`1/\nu`. It is in this
situation that the new option may be of great value.

In terms of partial differential equations, the typical problems for
which the stability limit detection option is appropriate are ODE
systems resulting from semi-discretized PDEs (i.e., PDEs discretized in
space) with advection and diffusion, but with advection dominating over
diffusion. Diffusion alone produces pure decay modes, while advection
tends to produce undamped oscillatory modes. A mix of the two with
advection dominant will have weakly damped oscillatory modes.

The STALD algorithm attempts to detect, in a direct manner, the
presence of a stability region boundary that is limiting the step sizes
in the presence of a weakly damped oscillation
:cite:p:`Hin:92`. The algorithm supplements (but differs
greatly from) the existing algorithms in CVODES for choosing step
size and order based on estimated local truncation errors. The STALD
algorithm works directly with history data that is readily available in
CVODES. If it concludes that the step size is in fact
stability-limited, it dictates a reduction in the method order,
regardless of the outcome of the error-based algorithm. The STALD
algorithm has been tested in combination with the VODE solver on
linear advection-dominated advection-diffusion problems
:cite:p:`Hin:95`, where it works well. The implementation in
CVODES has been successfully tested on linear and nonlinear
advection-diffusion problems, among others.

This stability limit detection option adds some computational overhead
to the CVODES solution. (In timing tests, these overhead costs have
ranged from 2% to 7% of the total, depending on the size and complexity
of the problem, with lower relative costs for larger problems.)
Therefore, it should be activated only when there is reasonable
expectation of modes in the user’s system for which it is appropriate.
In particular, if a CVODES solution with this option turned off
appears to take an inordinately large number of steps at orders 3-5 for
no apparent reason in terms of the solution time scale, then there is a
good chance that step sizes are being limited by stability, and that
turning on the option will improve the efficiency of the solution.

.. _CVODES.Mathematics.rootfinding:

Rootfinding
===========

The CVODES solver has been augmented to include a rootfinding
feature. This means that, while integrating the Initial Value Problem
:eq:`CVODES_ivp`, CVODES can also find the roots of a set of
user-defined functions :math:`g_i(t,y)` that depend both on :math:`t`
and on the solution vector :math:`y = y(t)`. The number of these root
functions is arbitrary, and if more than one :math:`g_i` is found to
have a root in any given interval, the various root locations are found
and reported in the order that they occur on the :math:`t` axis, in the
direction of integration.

Generally, this rootfinding feature finds only roots of odd
multiplicity, corresponding to changes in sign of :math:`g_i(t,y(t))`,
denoted :math:`g_i(t)` for short. If a user root function has a root of
even multiplicity (no sign change), it will probably be missed by
CVODES. If such a root is desired, the user should reformulate the
root function so that it changes sign at the desired root.

The basic scheme used is to check for sign changes of any :math:`g_i(t)`
over each time step taken, and then (when a sign change is found) to
hone in on the root(s) with a modified secant method
:cite:p:`HeSh:80`. In addition, each time :math:`g` is
computed, CVODES checks to see if :math:`g_i(t) = 0` exactly, and if
so it reports this as a root. However, if an exact zero of any
:math:`g_i` is found at a point :math:`t`, CVODES computes :math:`g`
at :math:`t + \delta` for a small increment :math:`\delta`, slightly
further in the direction of integration, and if any
:math:`g_i(t + \delta)=0` also, CVODES stops and reports an error.
This way, each time CVODES takes a time step, it is guaranteed that
the values of all :math:`g_i` are nonzero at some past value of
:math:`t`, beyond which a search for roots is to be done.

At any given time in the course of the time-stepping, after suitable
checking and adjusting has been done, CVODES has an interval
:math:`(t_{lo},t_{hi}]` in which roots of the :math:`g_i(t)` are to be
sought, such that :math:`t_{hi}` is further ahead in the direction of
integration, and all :math:`g_i(t_{lo}) \neq 0`. The endpoint
:math:`t_{hi}` is either :math:`t_n`, the end of the time step last
taken, or the next requested output time
:math:`t_{\text{out}}` if this comes sooner. The endpoint
:math:`t_{lo}` is either :math:`t_{n-1}`, the last output time
:math:`t_{\text{out}}` (if this occurred within the last
step), or the last root location (if a root was just located within this
step), possibly adjusted slightly toward :math:`t_n` if an exact zero
was found. The algorithm checks :math:`g_i` at :math:`t_{hi}` for zeros
and for sign changes in :math:`(t_{lo},t_{hi})`. If no sign changes were
found, then either a root is reported (if some :math:`g_i(t_{hi}) = 0`)
or we proceed to the next time interval (starting at :math:`t_{hi}`). If
one or more sign changes were found, then a loop is entered to locate
the root to within a rather tight tolerance, given by

.. math:: \tau = 100 * U * (|t_n| + |h|)~~~ (U = \mbox{unit roundoff}) ~.

Whenever sign changes are seen in two or more root functions, the one
deemed most likely to have its root occur first is the one with the
largest value of :math:`|g_i(t_{hi})|/|g_i(t_{hi}) - g_i(t_{lo})|`,
corresponding to the closest to :math:`t_{lo}` of the secant method
values. At each pass through the loop, a new value :math:`t_{mid}` is
set, strictly within the search interval, and the values of
:math:`g_i(t_{mid})` are checked. Then either :math:`t_{lo}` or
:math:`t_{hi}` is reset to :math:`t_{mid}` according to which
subinterval is found to include the sign change. If there is none in
:math:`(t_{lo},t_{mid})` but some :math:`g_i(t_{mid}) = 0`, then that
root is reported. The loop continues until
:math:`|t_{hi}-t_{lo}| < \tau`, and then the reported root location is
:math:`t_{hi}`.

In the loop to locate the root of :math:`g_i(t)`, the formula for
:math:`t_{mid}` is

.. math::

   t_{mid} = t_{hi} - (t_{hi} - t_{lo})
                g_i(t_{hi}) / [g_i(t_{hi}) - \alpha g_i(t_{lo})] ~,

where :math:`\alpha` is a weight parameter. On the first two passes
through the loop, :math:`\alpha` is set to :math:`1`, making
:math:`t_{mid}` the secant method value. Thereafter, :math:`\alpha` is
reset according to the side of the subinterval (low vs. high, i.e.,
toward :math:`t_{lo}` vs. toward :math:`t_{hi}`) in which the sign
change was found in the previous two passes. If the two sides were
opposite, :math:`\alpha` is set to 1. If the two sides were the same,
:math:`\alpha` is halved (if on the low side) or doubled (if on the high
side). The value of :math:`t_{mid}` is closer to :math:`t_{lo}` when
:math:`\alpha < 1` and closer to :math:`t_{hi}` when :math:`\alpha > 1`.
If the above value of :math:`t_{mid}` is within :math:`\tau/2` of
:math:`t_{lo}` or :math:`t_{hi}`, it is adjusted inward, such that its
fractional distance from the endpoint (relative to the interval size) is
between .1 and .5 (.5 being the midpoint), and the actual distance from
the endpoint is at least :math:`\tau/2`.


.. _CVODES.Mathematics.quad:

Pure Quadrature Integration
===========================

In many applications, and most notably during the backward integration phase of
an adjoint sensitivity analysis run (see :numref:`CVODES.Mathematics.ASA`)
it is of interest to compute integral quantities of the form

.. math::
   z(t) = \int_{t_0}^t q(\tau, y(\tau), p) \, \mathrm d\tau \, .
   :label: CVODES_QUAD

The most effective approach to compute :math:`z(t)` is to extend the original
problem with the additional ODEs (obtained by applying Leibnitz’s
differentiation rule):

.. math:: \dot z = q(t,y,p) \, , \quad z(t_0) = 0 \, .

Note that this is equivalent to using a quadrature method based on the
underlying linear multistep polynomial representation for :math:`y(t)`.

This can be done at the “user level” by simply exposing to CVODES the extended
ODE system :eq:`CVODES_ivp_p` + :eq:`CVODES_QUAD`. However, in the
context of an implicit integration solver, this approach is not desirable since
the nonlinear solver module will require the Jacobian (or Jacobian-vector
product) of this extended ODE. Moreover, since the additional states :math:`z`
do not enter the right-hand side of the ODE :eq:`CVODES_QUAD` and
therefore the right-hand side of the extended ODE system, it is much more
efficient to treat the ODE system :eq:`CVODES_QUAD` separately from the
original system :eq:`CVODES_ivp_p` by “taking out” the additional states
:math:`z` from the nonlinear system :eq:`CVODES_nonlinear` that must
be solved in the correction step of the LMM. Instead, “corrected” values
:math:`z^n` are computed explicitly as

.. math::

   z^n = - \frac{1}{\alpha_{n,0}} \left(
       h_n \beta_{n,0} q(t_n, y_n, p) + h_n \sum_{i=1}^{K_2} \beta_{n,i} \dot
       z^{n-i} + \sum_{i=1}^{K_1} \alpha_{n,i} z^{n-i} \right) \, ,

once the new approximation :math:`y^n` is available.

The quadrature variables :math:`z` can be optionally included in the error test,
in which case corresponding relative and absolute tolerances must be provided.


.. _CVODES.Mathematics.FSA:

Forward Sensitivity Analysis
============================

Typically, the governing equations of complex, large-scale models depend on
various parameters, through the right-hand side vector and/or through the vector
of initial conditions, as in :eq:`CVODES_ivp_p`. In addition to
numerically solving the ODEs, it may be desirable to determine the sensitivity
of the results with respect to the model parameters. Such sensitivity
information can be used to estimate which parameters are most influential in
affecting the behavior of the simulation or to evaluate optimization gradients
(in the setting of dynamic optimization, parameter estimation, optimal control,
etc.).

The *solution sensitivity* with respect to the model parameter :math:`p_i` is
defined as the vector :math:`s_i (t) = {\partial y(t)}/{\partial p_i}` and
satisfies the following *forward sensitivity equations* (or *sensitivity
equations* for short):

.. math::
   {{\dot s}_i}  = \frac{\partial f}{\partial y} s_i +
   \frac{\partial f}{\partial p_i} \, , \quad s_i(t_0)  = \frac{\partial
   y_{0}(p)}{\partial p_i} \, ,
   :label: CVODES_sens_eqns

obtained by applying the chain rule of differentiation to the original
ODEs :eq:`CVODES_ivp_p`.

When performing forward sensitivity analysis, CVODES carries out the time
integration of the combined system, :eq:`CVODES_ivp_p` and
:eq:`CVODES_sens_eqns`, by viewing it as an ODE system of size
:math:`N(N_s+1)`, where :math:`N_s` is the number of model parameters
:math:`p_i`, with respect to which sensitivities are desired (:math:`N_s \le
N_p`). However, major improvements in efficiency can be made by taking advantage
of the special form of the sensitivity equations as linearizations of the
original ODEs. In particular, for stiff systems, for which CVODES employs a
Newton iteration, the original ODE system and all sensitivity systems share the
same Jacobian matrix, and therefore the same iteration matrix :math:`M` in
:eq:`CVODES_Newtonmat`.

The sensitivity equations are solved with the same linear multistep formula that
was selected for the original ODEs and, if Newton iteration was selected, the
same linear solver is used in the correction phase for both state and
sensitivity variables. In addition, CVODES offers the option of including (*full
error control*) or excluding (*partial error control*) the sensitivity variables
from the local error test.


Forward sensitivity methods
---------------------------

In what follows we briefly describe three methods that have been proposed for
the solution of the combined ODE and sensitivity system for the vector
:math:`{\hat y} = [y, s_1, \ldots , s_{N_s}]`.

-  *Staggered Direct*

   In this approach :cite:p:`CaSt:85`, the nonlinear system :eq:`CVODES_nonlinear` is first solved and, once an acceptable numerical solution
   is obtained, the sensitivity variables at the new step are found by directly
   solving :eq:`CVODES_sens_eqns` after the (BDF or Adams)
   discretization is used to eliminate :math:`{\dot s}_i`. Although the system
   matrix of the above linear system is based on exactly the same information as
   the matrix :math:`M` in :eq:`CVODES_Newtonmat`, it must be
   updated and factored at every step of the integration, in contrast to an
   evalutaion of :math:`M` which is updated only occasionally. For problems with
   many parameters (relative to the problem size), the staggered direct method
   can outperform the methods described below :cite:p:`LPZ:99`. However, the
   computational cost associated with matrix updates and factorizations makes
   this method unattractive for problems with many more states than parameters
   (such as those arising from semidiscretization of PDEs) and is therefore not
   implemented in CVODES.

-  *Simultaneous Corrector*

   In this method :cite:p:`MaPe:97`, the discretization is applied
   simultaneously to both the original equations :eq:`CVODES_ivp_p` and
   the sensitivity systems :eq:`CVODES_sens_eqns` resulting in the
   following nonlinear system

   .. math::
      {\hat F}({\hat y}_n) \equiv
         {\hat y}_n - h_n\beta_{n,0} {\hat f}(t_n,\,{\hat y}_n) - {\hat a}_n = 0 \, ,

   where :math:`{\hat f} = [ f(t,y,p), \ldots, ({\partial f}/{\partial y})(t,y,p) s_i + ({\partial f}/{\partial p_i})(t,y,p), \ldots ]`,
   and :math:`{\hat a}_n` is comprised of the terms in the discretization that
   depend on the solution at previous integration steps. This combined nonlinear
   system can be solved using a modified Newton method as in :eq:`CVODES_Newton` by solving
   the corrector equation

   .. math::
      {\hat M}[{\hat y}_{n(m+1)}-{\hat y}_{n(m)}]=-{\hat F}({\hat y}_{n(m)})
      :label: CVODES_Newton_sim

   at each iteration, where

   .. math::

      {\hat M} =
          \begin{bmatrix}
            M                &        &        &        &   \\ - \gamma J_1
            & M      &        &        &   \\ - \gamma J_2     & 0      & M
            &        &   \\
              \vdots         & \vdots & \ddots & \ddots &   \\
            - \gamma J_{N_s} & 0      & \ldots & 0      & M
          \end{bmatrix} \, ,

   :math:`M` is defined as in :eq:`CVODES_Newtonmat`, and
   :math:`J_i = \dfrac{\partial}{\partial y}\left[ \left(\dfrac{\partial f}{\partial y}\right) s_i + \left(\dfrac{\partial f}{\partial p_i}\right) \right]`. It
   can be shown that 2-step quadratic convergence can be retained by using only
   the block-diagonal portion of :math:`{\hat M}` in the corrector equation
   :eq:`CVODES_Newton_sim`. This results in a decoupling that
   allows the reuse of :math:`M` without additional matrix factorizations.
   However, the products :math:`\left(\dfrac{\partial f}{\partial y}\right)s_i` and the vectors
   :math:`\dfrac{\partial f}{\partial p_i}`
   must still be reevaluated at each step of the iterative process
   :eq:`CVODES_Newton_sim` to update the sensitivity portions of
   the residual :math:`{\hat G}`.

-  *Staggered corrector*

   In this approach :cite:p:`FTB:97`, as in the staggered direct method, the
   nonlinear system :eq:`CVODES_nonlinear` is solved first using the
   Newton iteration :eq:`CVODES_Newton`. Then a separate Newton
   iteration is used to solve the sensitivity system :eq:`CVODES_sens_eqns`:

   .. math::
      \begin{gathered}
          M [s_{i}^{n(m+1)} - s_{i}^{n(m)}]= \\ - \left[ s_{i}^{n(m)} - \gamma
          \left( \dfrac{\partial f}{\partial y} (t_n , y^n, p) s_{i}^{n(m)} + \dfrac{\partial f}{\partial p_i} (t_n , y^n , p)
          \right) -a_{i,n} \right] \, ,
      \end{gathered}
      :label: CVODES_stgr_iterations

   where :math:`a_{i,n}=\sum_{j>0}(\alpha_{n,j}s_{i}^{n-j}+h_n\beta_{n,j}{{\dot s}_i}^{n-j})`.
   In other words, a modified Newton iteration is used to solve a
   linear system. In this approach, the vectors :math:`({\partial f}/{\partial p_i})` need be updated
   only once per integration step, after the state correction phase :eq:`CVODES_Newton` has converged.
   Note also that Jacobian-related data can be reused at all iterations :eq:`CVODES_stgr_iterations` to
   evaluate the products :math:`({\partial f}/{\partial y}) s_i`.

CVODES implements the simultaneous corrector method and two flavors of the
staggered corrector method which differ only if the sensitivity variables are
included in the error control test. In the *full error control* case, the first
variant of the staggered corrector method requires the convergence of the
iterations :eq:`CVODES_stgr_iterations` for all :math:`N_s`
sensitivity systems and then performs the error test on the sensitivity
variables. The second variant of the method will perform the error test for each
sensitivity vector :math:`s_i, (i=1,2,\ldots,N_s`) individually, as they pass
the convergence test. Differences in performance between the two variants may
therefore be noticed whenever one of the sensitivity vectors :math:`s_i` fails a
convergence or error test.

An important observation is that the staggered corrector method, combined with a
Krylov linear solver, effectively results in a staggered direct method. Indeed,
the Krylov solver requires only the action of the matrix :math:`M` on a vector
and this can be provided with the current Jacobian information. Therefore, the
modified Newton procedure :eq:`CVODES_stgr_iterations` will
theoretically converge after one iteration.

Selection of the absolute tolerances for sensitivity variables
--------------------------------------------------------------

If the sensitivities are included in the error test, CVODES provides an
automated estimation of absolute tolerances for the sensitivity variables based
on the absolute tolerance for the corresponding state variable. The relative
tolerance for sensitivity variables is set to be the same as for the state
variables. The selection of absolute tolerances for the sensitivity variables is
based on the observation that the sensitivity vector :math:`s_i` will have units
of :math:`[y]/[p_i]`. With this, the absolute tolerance for the :math:`j`-th
component of the sensitivity vector :math:`s_i` is set to
:math:`{\mbox{atol}_j}/{|{\bar p}_i|}`, where :math:`\mbox{atol}_j` are the absolute
tolerances for the state variables and :math:`\bar p` is a vector of scaling
factors that are dimensionally consistent with the model parameters :math:`p`
and give an indication of their order of magnitude. This choice of relative and
absolute tolerances is equivalent to requiring that the weighted
root-mean-square norm of the sensitivity vector :math:`s_i` with weights based
on :math:`s_i` be the same as the weighted root-mean-square norm of the vector
of scaled sensitivities :math:`{\bar s}_i = |{\bar p}_i| s_i` with weights based
on the state variables (the scaled sensitivities :math:`{\bar s}_i` being
dimensionally consistent with the state variables). However, this choice of
tolerances for the :math:`s_i` may be a poor one, and the user of CVODES can
provide different values as an option.


Evaluation of the sensitivity right-hand side
---------------------------------------------

There are several methods for evaluating the right-hand side of the sensitivity
systems :eq:`CVODES_sens_eqns`: analytic evaluation, automatic
differentiation, complex-step approximation, and finite differences (or
directional derivatives). CVODES provides all the software hooks for
implementing interfaces to automatic differentiation (AD) or complex-step
approximation; future versions will include a generic interface to AD-generated
functions. At the present time, besides the option for analytical sensitivity
right-hand sides (user-provided), CVODES can evaluate these quantities using
various finite difference-based approximations to evaluate the terms
:math:`({\partial f}/{\partial y}) s_i` and :math:`({\partial f}/{\partial p_i})`, or using directional derivatives to
evaluate :math:`\left[ ({\partial f}/{\partial y}) s_i + ({\partial f}/{\partial p_i}) \right]`. As is typical for
finite differences, the proper choice of perturbations is a delicate matter.
CVODES takes into account several problem-related features: the relative ODE
error tolerance :math:`\mbox{rtol}`, the machine unit roundoff :math:`U`,
the scale factor :math:`{\bar p}_i`, and the weighted root-mean-square norm of
the sensitivity vector :math:`s_i`.

Using central finite differences as an example, the two terms :math:`({\partial
f}/{\partial y}) s_i` and :math:`{\partial f}/{\partial p_i}` in the right-hand
side of :eq:`CVODES_sens_eqns` can be evaluated either separately:

.. math::
   \frac{\partial f}{\partial y} s_i \approx \frac{f(t, y+\sigma_y s_i, p)-
      f(t, y-\sigma_y s_i, p)}{2\,\sigma_y} \, , \\
   :label: CVODES_fd2_1

.. math::
   \frac{\partial f}{\partial p_i} \approx \frac{f(t,y,p + \sigma_i e_i)-
      f(t,y,p - \sigma_i e_i)}{2\,\sigma_i} \, , \\
   :label: CVODES_fd2_2

.. math::
   \sigma_i = |{\bar p}_i| \sqrt{\max( \mbox{rtol} , U)} \, , \quad
   \sigma_y = \frac{1}{\max(1/\sigma_i, \|s_i\|/|{\bar p}_i|)} \, ,

or simultaneously:

.. math::
   \begin{gathered}
     \frac{\partial f}{\partial y} s_i + \frac{\partial f}{\partial p_i} \approx
     \frac{f(t, y+\sigma s_i, p + \sigma e_i) -
       f(t, y-\sigma s_i, p - \sigma e_i)}{2\,\sigma} \, , \\
     \sigma = \min(\sigma_i, \sigma_y) \, , \nonumber\end{gathered}
   :label: CVODES_dd2

or by adaptively switching between :eq:`CVODES_fd2_1` + :eq:`CVODES_fd2_2` and :eq:`CVODES_dd2`, depending on the relative size of the finite
difference increments :math:`\sigma_i` and :math:`\sigma_y`. In the adaptive
scheme, if :math:`\rho = \max(\sigma_i/\sigma_y,\sigma_y/\sigma_i)`, we use
separate evaluations if :math:`\rho > \rho_{max}` (an input value), and
simultaneous evaluations otherwise.

These procedures for choosing the perturbations (:math:`\sigma_i`,
:math:`\sigma_y`, :math:`\sigma`) and switching between finite difference and
directional derivative formulas have also been implemented for one-sided
difference formulas. Forward finite differences can be applied to
:math:`({\partial f}/{\partial y}) s_i` and :math:`{\partial f}/{\partial p_i}`
separately, or the single directional derivative formula

.. math::

   \dfrac{\partial f}{\partial y} s_i + \dfrac{\partial f}{\partial p_i} \approx \frac{f(t, y+\sigma s_i, p + \sigma e_i) - f(t, y,
   p)}\sigma

can be used. In CVODES, the default value of :math:`\rho_{max}=0` indicates the use
of the second-order centered directional derivative formula :eq:`CVODES_dd2` exclusively.
Otherwise, the magnitude of :math:`\rho_{max}` and its
sign (positive or negative) indicates whether this switching is done with regard
to (centered or forward) finite differences, respectively.


Quadratures depending on forward sensitivities
----------------------------------------------

If pure quadrature variables are also included in the problem definition (see
:numref:`CVODES.Mathematics.quad`), CVODES does *not* carry their sensitivities
automatically. Instead, we provide a more general feature through which
integrals depending on both the states :math:`y` of :eq:`CVODES_ivp_p`
and the state sensitivities :math:`s_i` of :eq:`CVODES_sens_eqns`
can be evaluated. In other words, CVODES provides support for computing
integrals of the form:

.. math:: \bar z(t) = \int_{t_0}^t \bar q(\tau, y(\tau), s_1(\tau), \ldots, s_{N_p}(\tau),p) \, \mathrm d\tau \, .

If the sensitivities of the quadrature variables :math:`z` of :eq:`CVODES_QUAD` are desired, these can then be computed by using:

.. math:: \bar q_i = q_y s_i + q_{p_i} \, , \quad i = 1,\ldots,N_p \, ,

as integrands for :math:`\bar{z}`, where :math:`q_y` and :math:`q_p` are the
partial derivatives of the integrand function :math:`q` of :eq:`CVODES_QUAD`.

As with the quadrature variables :math:`z`, the new variables :math:`\bar z` are
also excluded from any nonlinear solver phase and “corrected” values :math:`\bar
z^n` are obtained through explicit formulas.


.. _CVODES.Mathematics.ASA:

Adjoint Sensitivity Analysis
============================

In the *forward sensitivity approach* described in the previous section,
obtaining sensitivities with respect to :math:`N_s` parameters is roughly
equivalent to solving an ODE system of size :math:`(1+N_s) N`. This can become
prohibitively expensive, especially for large-scale problems, if sensitivities
with respect to many parameters are desired. In this situation, the *adjoint
sensitivity method* is a very attractive alternative, provided that we do not
need the solution sensitivities :math:`s_i`, but rather the gradients with
respect to model parameters of a relatively few derived functionals of the
solution. In other words, if :math:`y(t)` is the solution of :eq:`CVODES_ivp_p`, we
wish to evaluate the gradient :math:`{\mathrm dG}/{\mathrm dp}` of

.. math::
   G(p) = \int_{t_0}^T g(t, y, p) \mathrm dt \, ,
   :label: CVODES_G

or, alternatively, the gradient :math:`{\mathrm dg}/{\mathrm dp}` of the function
:math:`g(t, y, p)` at the final time :math:`T`. The function :math:`g` must be smooth enough
that :math:`\partial g / \partial y` and :math:`\partial g / \partial p` exist
and are bounded.

In what follows, we only sketch the analysis for the sensitivity problem for
both :math:`G` and :math:`g`. For details on the derivation see
:cite:p:`CLPS:03`. Introducing a Lagrange multiplier :math:`\lambda`, we form
the augmented objective function

.. math::

   I(p) = G(p) - \int_{t_0}^T \lambda^*
   \left( {\dot y} - f(t,y,p)\right) \mathrm dt \, ,

where :math:`*` denotes the conjugate transpose. The gradient of :math:`G` with respect to :math:`p` is

.. math::

   \frac{\mathrm dG}{\mathrm dp} = \frac{\mathrm dI}{\mathrm dp}
   =\int_{t_0}^T(g_p + g_y s) \mathrm dt - \int_{t_0}^T
   \lambda^* \left( {\dot s} - f_y s - f_p \right)\mathrm dt \, ,

where subscripts on functions :math:`f` or :math:`g` are used to denote partial
derivatives and :math:`s = [s_1,\ldots,s_{N_s}]` is the matrix of solution
sensitivities. Applying integration by parts to the term :math:`\lambda^* {\dot
s}`, and by requiring that :math:`\lambda` satisfy

.. math::
   \begin{split}
   &{\dot\lambda} = -\left( \frac{\partial f}{\partial y} \right)^* \lambda -
   \left( \frac{\partial g}{\partial y} \right)^* \\
   &\lambda(T) = 0 \, ,
   \end{split}
   :label: CVODES_adj_eqns

the gradient of :math:`G` with respect to :math:`p` is nothing but

.. math::
   \frac{\mathrm dG}{\mathrm dp} = \lambda^*(t_0) s(t_0) +
   \int_{t_0}^T \left( g_p + \lambda^* f_p \right) \mathrm dt \, .
   :label: CVODES_dgdp_1

The gradient of :math:`g(T,y,p)` with respect to :math:`p` can be then obtained
by using the Leibnitz differentiation rule. Indeed, from :eq:`CVODES_G`,

.. math:: \frac{\mathrm dg}{\mathrm dp}(T) = \frac{\mathrm d}{\mathrm dT}\frac{\mathrm dG}{\mathrm dp}

and therefore, taking into account that :math:`dG/dp` in :eq:`CVODES_dgdp_1` depends on
:math:`T` both through the upper integration limit and through :math:`\lambda`,
and that :math:`\lambda(T) = 0`,

.. math::
   \frac{\mathrm dg}{\mathrm dp}(T) = \mu^*(t_0) s(t_0) + g_p(T) +
   \int_{t_0}^T \mu^* f_p \mathrm dt \, ,
   :label: CVODES_dgdp_2

where :math:`\mu` is the sensitivity of :math:`\lambda` with respect to the
final integration limit :math:`T`. Thus :math:`\mu` satisfies the following
equation, obtained by taking the total derivative with respect to :math:`T` of
:eq:`CVODES_adj_eqns`:

.. math::
   \begin{split}
   &{\dot\mu} = -\left( \frac{\partial f}{\partial y} \right)^* \mu \\
   &\mu(T) = \left( \frac{\partial g}{\partial y} \right)^*_{t=T} \, .
   \end{split}
   :label: CVODES_adj1_eqns

The final condition on :math:`\mu(T)` follows from
:math:`(\partial\lambda/\partial t) + (\partial\lambda/\partial T) = 0` at
:math:`T`, and therefore, :math:`\mu(T) = -{\dot\lambda}(T)`.

The first thing to notice about the adjoint system :eq:`CVODES_adj_eqns` is that there
is no explicit specification of the parameters :math:`p`; this implies that,
once the solution :math:`\lambda` is found, the formula :eq:`CVODES_dgdp_1` can then be
used to find the gradient of :math:`G` with respect to any of the parameters
:math:`p`. The same holds true for the system :eq:`CVODES_adj1_eqns` and the formula
:eq:`CVODES_dgdp_2` for gradients of :math:`g(T,y,p)`. The second important remark is
that the adjoint systems :eq:`CVODES_adj_eqns` and :eq:`CVODES_adj1_eqns` are terminal value
problems which depend on the solution :math:`y(t)` of the original IVP
:eq:`CVODES_ivp_p`. Therefore, a procedure is needed for providing the states :math:`y`
obtained during a forward integration phase of :eq:`CVODES_ivp_p` to CVODES during the
backward integration phase of :eq:`CVODES_adj_eqns` or :eq:`CVODES_adj1_eqns`. The approach
adopted in CVODES, based on *checkpointing*, is described below.

.. _CVODES.Mathematics.Checkpointing:

Checkpointing scheme
====================

During the backward integration, the evaluation of the right-hand side of the
adjoint system requires, at the current time, the states :math:`y` which were
computed during the forward integration phase. Since CVODES implements
variable-step integration formulas, it is unlikely that the states will be
available at the desired time and so some form of interpolation is needed. The
CVODES implementation being also variable-order, it is possible that during the
forward integration phase the order may be reduced as low as first order, which
means that there may be points in time where only :math:`y` and :math:`{\dot y}`
are available. These requirements therefore limit the choices for possible
interpolation schemes. CVODES implements two interpolation methods: a cubic
Hermite interpolation algorithm and a variable-degree polynomial interpolation
method which attempts to mimic the BDF interpolant for the forward integration.

However, especially for large-scale problems and long integration intervals, the
number and size of the vectors :math:`y` and :math:`{\dot y}` that would need to
be stored make this approach computationally intractable. Thus, CVODES settles
for a compromise between storage space and execution time by implementing a
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
integrated from :math:`T` to :math:`t_0` going from one checkpoint to the
previous one. The backward integration from checkpoint :math:`i+1` to checkpoint
:math:`i` is preceded by a forward integration from :math:`i` to :math:`i+1`
during which the :math:`N_d` vectors :math:`y` (and, if necessary :math:`{\dot
y}`) are generated and stored in memory for interpolation
(see :numref:`CVODES.Mathematics.Checkpointing.Figure`).

.. note::
   The degree of the interpolation polynomial is always that of the current BDF
   order for the forward interpolation at the first point to the right of the
   time at which the interpolated value is sought (unless too close to the
   :math:`i`-th checkpoint, in which case it uses the BDF order at the
   right-most relevant point). However, because of the FLC BDF implementation
   :numref:`CVODES.Mathematics.ivp_sol`, the resulting interpolation
   polynomial is only an approximation to the underlying BDF interpolant.

   The Hermite cubic interpolation option is present because it was implemented
   chronologically first and it is also used by other adjoint solvers (e.g.
   DASPKADJOINT. The variable-degree polynomial is more memory-efficient (it
   requires only half of the memory storage of the cubic Hermite interpolation)
   and is more accurate. The accuracy differences are minor when using BDF
   (since the maximum method order cannot exceed 5), but can be significant for
   the Adams method for which the order can reach 12.


.. _CVODES.Mathematics.Checkpointing.Figure:

.. figure:: /figs/cvodes/ckpnt.png
   :align: center
   :alt: Illustration of the checkpointing algorithm for generation of the forward solution during the integration of the adjoint system.

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
the solution of :eq:`CVODES_ivp_p`), the total cost of an adjoint sensitivity
computation can be as low as one forward plus one backward integration. In
addition, CVODES provides the capability of reusing a set of checkpoints for
multiple backward integrations, thus allowing for efficient computation of
gradients of several functionals :eq:`CVODES_G`.

Finally, we note that the adjoint sensitivity module in CVODES provides the
necessary infrastructure to integrate backwards in time any ODE terminal value
problem dependent on the solution of the IVP :eq:`CVODES_ivp_p`, including adjoint
systems :eq:`CVODES_adj_eqns` or :eq:`CVODES_adj1_eqns`, as well as any other quadrature ODEs
that may be needed in evaluating the integrals in :eq:`CVODES_dgdp_1` or :eq:`CVODES_dgdp_2`. In
particular, for ODE systems arising from semi-discretization of time-dependent
PDEs, this feature allows for integration of either the discretized adjoint PDE
system or the adjoint of the discretized PDE.

.. _CVODES.Mathematics.hess_sensi:

Second-order sensitivity analysis
=================================

In some applications (e.g., dynamically-constrained optimization) it may be
desirable to compute second-order derivative information. Considering the ODE
problem :eq:`CVODES_ivp_p` and some model output functional, :math:`g(y)` then the
Hessian :math:`d^2g/dp^2` can be obtained in a forward sensitivity analysis
setting as

.. math:: \frac{\mathrm d^2 g}{\mathrm d p^2} = \left(g_y \otimes I_{N_p} \right ) y_{pp} + y_p^T g_{yy} y_p \, ,

where :math:`\otimes` is the Kronecker product. The second-order sensitivities
are solution of the matrix ODE system:

.. math::

   \begin{split}
       & {\dot y}_{pp} = \left( f_y \otimes I_{N_p} \right) \cdot y_{pp} +
       \left( I_N \otimes y_p^T \right) \cdot f_{yy} y_p \\
       & y_{pp}(t_0) = \frac{\partial^2 y_0}{\partial p^2} \, ,
     \end{split}

where :math:`y_p` is the first-order sensitivity matrix, the solution of
:math:`N_p` systems :eq:`CVODES_sens_eqns`, and :math:`y_{pp}` is a third-order tensor.
It is easy to see that, except for situations in which the number of parameters
:math:`N_p` is very small, the computational cost of this so-called
*forward-over-forward* approach is exorbitant as it requires the solution of
:math:`N_p + N_p^2` additional ODE systems of the same dimension :math:`N` as
:eq:`CVODES_ivp_p`.

.. note::
   For the sake of simplifity in presentation, we do not include explicit
   dependencies of :math:`g` on time :math:`t` or parameters :math:`p`.
   Moreover, we only consider the case in which the dependency of the original
   ODE :eq:`CVODES_ivp_p` on the parameters :math:`p` is through its initial conditions
   only. For details on the derivation in the general case,
   see :cite:p:`OzBa:05`.

A much more efficient alternative is to compute Hessian-vector products using a
so-called *forward-over-adjoint* approach. This method is based on using the
same “trick” as the one used in computing gradients of pointwise functionals
with the adjoint method, namely applying a formal directional forward derivation
to one of the gradients of :eq:`CVODES_dgdp_1` or :eq:`CVODES_dgdp_2`. With that, the cost of
computing a full Hessian is roughly equivalent to the cost of computing the
gradient with forward sensitivity analysis. However, Hessian-vector products can
be cheaply computed with one additional adjoint solve. Consider for example,
:math:`G(p) = \int_{t_0}^{t_f} g(t,y) \, \mathrm dt`. It can be shown that the product
between the Hessian of :math:`G` (with respect to the parameters :math:`p`) and
some vector :math:`u` can be computed as

.. math::

   \frac{\partial^2 G}{\partial p^2} u =
     \left[ \left(\lambda^T \otimes I_{N_p} \right) y_{pp}u + y_p^T \mu \right]_{t=t_0} \, ,

where :math:`\lambda`, :math:`\mu`, and :math:`s` are solutions of

.. math::

   \begin{split}
       &-\dot\mu = f_y^T\mu + \left(\lambda^T \otimes I_n \right) f_{yy} s + g_{yy} s\, ; \quad \mu(t_f) = 0 \\
       &-\dot\lambda = f_y^T\lambda + g_y^T \, ; \quad \lambda(t_f) = 0 \\
       &\dot s = f_y s\, ; \quad s(t_0) = y_{0p} u
     \end{split}

In the above equation, :math:`s = y_p u` is a linear combination of the columns
of the sensitivity matrix :math:`y_p`. The *forward-over-adjoint* approach
hinges crucially on the fact that :math:`s` can be computed at the cost of a
forward sensitivity analysis with respect to a single parameter (the last ODE
problem above) which is possible due to the linearity of the forward sensitivity
equations :eq:`CVODES_sens_eqns`.

Therefore, the cost of computing the Hessian-vector product is roughly that of
two forward and two backward integrations of a system of ODEs of size :math:`N`.
For more details, including the corresponding formulas for a pointwise model
functional output, see :cite:p:`OzBa:05`.

To allow the *foward-over-adjoint* approach described above, CVODES provides
support for:

-  the integration of multiple backward problems depending on the same
   underlying forward problem :eq:`CVODES_ivp_p`, and

-  the integration of backward problems and computation of backward
   quadratures depending on both the states :math:`y` and forward
   sensitivities (for this particular application, :math:`s`) of the original
   problem :eq:`CVODES_ivp_p`.
