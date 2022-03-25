.. ----------------------------------------------------------------
   SUNDIALS Copyright Start
   Copyright (c) 2002-2022, Lawrence Livermore National Security
   and Southern Methodist University.
   All rights reserved.

   See the top-level LICENSE and NOTICE files for details.

   SPDX-License-Identifier: BSD-3-Clause
   SUNDIALS Copyright End
   ----------------------------------------------------------------

.. _CVODE.Mathematics:

***************************
Mathematical Considerations
***************************

CVODE solves ODE initial value problems (IVPs) in real
:math:`N`-space, which we write in the abstract form

.. math::
   \dot{y} = f(t,y) \, ,\quad y(t_0) = y_0 \,
   :label: CVODE_ivp

where :math:`y \in \mathbb{R}^N` and
:math:`f: \mathbb{R} \times \mathbb{R}^N \rightarrow \mathbb{R}^N`.
Here we use :math:`\dot{y}` to denote :math:`\mathrm dy/\mathrm dt`. While we use
:math:`t` to denote the independent variable, and usually this is time,
it certainly need not be. CVODE solves both stiff and nonstiff
systems. Roughly speaking, stiffness is characterized by the presence of
at least one rapidly damped mode, whose time constant is small compared
to the time scale of the solution itself. Additionally, for problems :eq:`CVODE_ivp`
where the analytical solution :math:`y(t)` satisfies an implicit constraint
:math:`g(t,y)=0` (including the initial condition, :math:`g(t_0,y_0)=0`) for
:math:`g(t,y): \mathbb{R} \times \mathbb{R}^N \rightarrow \mathbb{R}^{M}` with
:math:`M<N`,  CVODE may be configured to explicitly enforce these constraints
via solving the modified problem

.. math::
   \begin{aligned}
      \dot{y} &= f(t,y) \, ,\quad y(t_0) = y_0 \, ,\\
      0 &= g(t,y).
   \end{aligned}
   :label: CVODE_ivp_constr


.. _CVODE.Mathematics.ivp_sol:

IVP solution
============

The methods used in CVODE are variable-order, variable-step
multistep methods, based on formulas of the form

.. math::
   \sum_{i = 0}^{K_1} \alpha_{n,i} y^{n-i} +
        h_n \sum_{i = 0}^{K_2} \beta_{n,i} {\dot{y}}^{n-i} = 0 \, .
   :label: CVODE_lmm

Here the :math:`y^n` are computed approximations to :math:`y(t_n)`, and
:math:`h_n = t_n - t_{n-1}` is the step size. The user of CVODE must
choose appropriately one of two multistep methods. For nonstiff
problems, CVODE includes the Adams-Moulton formulas, characterized
by :math:`K_1 = 1` and :math:`K_2 = q-1` above, where the order
:math:`q` varies between :math:`1` and :math:`12`. For stiff problems,
CVODE includes the Backward Differentiation Formulas (BDF) in
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
   :label: CVODE_nonlinear

or as a fixed-point problem

.. math::
   G(y^n) \equiv h_n \beta_{n,0} f(t_n,y^n) + a_n = y^n \, .
   :label: CVODE_nonlinear_fixedpoint

where
:math:`a_n\equiv\sum_{i>0}(\alpha_{n,i}y^{n-i}+h_n\beta_{n,i} {\dot{y}}^{n-i})`.
CVODE provides several nonlinear solver choices as well as the
option of using a user-defined nonlinear solver (see
:numref:`SUNNonlinSol`). By default CVODE solves :eq:`CVODE_nonlinear` with a
*Newton iteration* which requires the solution of linear systems

.. math::
   M [y^{n(m+1)} - y^{n(m)}] = -F(y^{n(m)}) \, ,
   :label: CVODE_Newton

in which

.. math::
   M \approx I - \gamma J \, ,
   \quad J = \partial f / \partial y \, ,
   \quad \mbox{and} \quad
   \gamma = h_n \beta_{n,0} \, .
   :label: CVODE_Newtonmat

The exact variation of the Newton iteration depends on the choice of linear
solver and is discussed below and in :numref:`SUNNonlinSol.Newton`. For nonstiff
systems, a *fixed-point iteration* (previously referred to as a functional
iteration in this guide) solving :eq:`CVODE_nonlinear_fixedpoint` is also
available. This involves evaluations of :math:`f` only and can (optionally) use
Anderson’s method :cite:p:`Anderson65, Walker-Ni09, Fang-Saad09, LWWY11` to
accelerate convergence (see :numref:`SUNNonlinSol.FixedPoint` for more details).
For any nonlinear solver, the initial guess for the iteration is a predicted
value :math:`y^{n(0)}` computed explicitly from the available history data.

For nonlinear solvers that require the solution of the linear system
:eq:`CVODE_Newton` (e.g., the default Newton iteration),
CVODE provides several linear solver choices, including the option
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

In addition, CVODE also provides a linear solver module which only
uses a diagonal approximation of the Jacobian matrix.

..
   Note that the dense, band, and sparse direct linear solvers can only be
   used with the serial and threaded vector representations. The diagonal
   solver can be used with any vector representation.

In the process of controlling errors at various levels, CVODE uses a
weighted root-mean-square norm, denoted
:math:`|\cdot|_{\text{WRMS}}`, for all error-like
quantities. The multiplicative weights used are based on the current
solution and on the relative and absolute tolerances input by the user,
namely

.. math::
   W_i = 1 / [\text{rtol} \cdot |y_i| + \text{atol}_i ] \, .
   :label: CVODE_errwt

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
:math:`y^n` denote the exact solution of :eq:`CVODE_nonlinear`, we want to ensure that the iteration
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
:eq:`CVODE_errwt`. In the dense case, this scheme requires
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
   :label: CVODE_jacobv

The increment :math:`\sigma` is :math:`1/\|v\|`, so that
:math:`\sigma v` has norm 1.

A critical part of CVODE — making it an ODE “solver” rather than
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
failures, CVODE returns to the user with a give-up message.

In addition to adjusting the step size to meet the local error test,
CVODE periodically adjusts the order, with the goal of maximizing
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

The various algorithmic features of CVODE described above, as
inherited from VODE and VODPK, are documented in
:cite:p:`BBH:89,Byr:92,Hin:00`. They are also summarized in
:cite:p:`HBGLSSW:05`.

CVODE permits the user to impose optional inequality constraints on
individual components of the solution vector :math:`y`. Any of the
following four constraints can be imposed: :math:`y_i > 0`,
:math:`y_i < 0`, :math:`y_i \geq 0`, or :math:`y_i \leq 0`. The
constraint satisfaction is tested after a successful nonlinear system
solution. If any constraint fails, we declare a convergence failure of
the Newton iteration and reduce the step size. Rather than cutting the
step size by some arbitrary factor, CVODE estimates a new step size
:math:`h'` using a linear approximation of the components in :math:`y`
that failed the constraint test (including a safety factor of
:math:`0.9` to cover the strict inequality case). If a step fails to
satisfy the constraints repeatedly within a step attempt or fails with
the minimum step size then the integration is halted and an error is
returned. In this case the user may need to employ other strategies as
discussed in :numref:`CVODE.Usage.CC.callable_fct_sim.cvtolerances` to satisfy
the inequality constraints.

Normally, CVODE takes steps until a user-defined output value
:math:`t = t_{\text{out}}` is overtaken, and then it
computes :math:`y(t_{\text{out}})` by interpolation.
However, a “one step” mode option is available, where control returns to
the calling program after each step. There are also options to force
CVODE not to integrate past a given stopping point
:math:`t = t_{\text{stop}}`.

.. _CVODE.Mathematics.constraints:

IVPs with constraints
=====================

For IVPs whose analytical solutions implicitly satisfy constraints as
in :eq:`CVODE_ivp_constr`, CVODE ensures that the solution satisfies
the constraint equation by projecting a successfully computed time step
onto the invariant manifold. As discussed in
:cite:p:`eich1993convergence` and
:cite:p:`shampine1999conservation`, this approach reduces the
error in the solution and retains the order of convergence of the
numerical method. Therefore, in an attempt to advance the solution to a
new point in time (i.e., taking a new integration step), CVODE
performs the following operations:

#. predict solution

#. solve nonlinear system and correct solution

#. project solution

#. test error

#. select order and step size for next step

and includes several recovery attempts in case there are convergence
failures (or difficulties) in the nonlinear solver or in the projection
step, or if the solution fails to satisfy the error test. Note that at
this time projection is only supported with BDF methods and the
projection function must be user-defined. See :numref:`CVODE.Usage.CC.cvprojinit` and
:c:func:`CVodeSetProjFn` for more information on providing a
projection function to CVODE.

When using a coordinate projection method the solution :math:`y_n` is
obtained by projecting (orthogonally or otherwise) the solution
:math:`\tilde{y}_n` from step 2 above onto
the manifold given by the constraint. As such :math:`y_n` is computed as
the solution of the nonlinear constrained least squares problem

.. math::
   \begin{split}
     \text{minimize}   &\quad \| y_n - \tilde{y}_n \| \\
     \text{subject to} &\quad g(t_n,y_n) = 0.
   \end{split}
   :label: CVODE_proj

The solution of :eq:`CVODE_proj` can be computed iteratively with
a Newton method. Given an initial guess :math:`y_n^{(0)}` the iterations
are computed as

.. math:: y_n^{(i+1)} = y_n^{(i)} + \delta y_n^{(i)}

where the increment :math:`\delta y_n^{(i)}` is the solution of the
least-norm problem

.. math::
   \begin{split}
       \text{minimize}   &\quad \| \delta y_n^{(i)} \| \\
       \text{subject to} &\quad G(t_n,y_n^{(i)}) \; \delta y_n^{(i)} = -g(t_n,y_n^{(i)})
   \end{split}
   :label: CVODE_leastnorm

where :math:`G(t,y) = \partial g(t,y) / \partial y`.

If the projected solution satisfies the error test then the step is
accepted and the correction to the unprojected solution,
:math:`\Delta_p = y_n - \tilde{y}_n`, is used to update the Nordsieck
history array for the next step.

.. _CVODE.Mathematics.preconditioning:

Preconditioning
===============

When using a nonlinear solver that requires the solution of the linear
system (:numref:`SUNNonlinSol.Newton`) (e.g., the default Newton
iteration), CVODE makes repeated use of a linear solver to solve
linear systems of the form :math:`M x = - r`, where :math:`x` is a
correction vector and :math:`r` is a residual vector. If this linear
system solve is done with one of the scaled preconditioned iterative
linear solvers supplied with SUNDIALS, these solvers are rarely
successful if used without preconditioning; it is generally necessary to
precondition the system in order to obtain acceptable efficiency. A
system :math:`A x = b` can be preconditioned on the left, as
:math:`(P^{-1}A) x = P^{-1} b`; on the right, as
:math:`(A P^{-1}) P x = b`; or on both sides, as
:math:`(P_L^{-1} A P_R^{-1}) P_R x = P_L^{-1}b`. The Krylov method is
then applied to a system with the matrix :math:`P^{-1}A`, or
:math:`AP^{-1}`, or :math:`P_L^{-1} A P_R^{-1}`, instead of :math:`A`.
In order to improve the convergence of the Krylov iteration, the
preconditioner matrix :math:`P`, or the product :math:`P_L P_R` in the
last case, should in some sense approximate the system matrix :math:`A`.
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

Typical preconditioners used with CVODE are based on approximations
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

.. _CVODE.Mathematics.stablimit:

BDF stability limit detection
=============================

CVODE includes an algorithm, STALD (STAbility Limit Detection),
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
CVODE. If it concludes that the step size is in fact
stability-limited, it dictates a reduction in the method order,
regardless of the outcome of the error-based algorithm. The STALD
algorithm has been tested in combination with the VODE solver on
linear advection-dominated advection-diffusion problems
:cite:p:`Hin:95`, where it works well. The implementation in
CVODE has been successfully tested on linear and nonlinear
advection-diffusion problems, among others.

This stability limit detection option adds some computational overhead
to the CVODES solution. (In timing tests, these overhead costs have
ranged from 2% to 7% of the total, depending on the size and complexity
of the problem, with lower relative costs for larger problems.)
Therefore, it should be activated only when there is reasonable
expectation of modes in the user’s system for which it is appropriate.
In particular, if a CVODE solution with this option turned off
appears to take an inordinately large number of steps at orders 3-5 for
no apparent reason in terms of the solution time scale, then there is a
good chance that step sizes are being limited by stability, and that
turning on the option will improve the efficiency of the solution.

.. _CVODE.Mathematics.rootfinding:

Rootfinding
===========

The CVODE solver has been augmented to include a rootfinding
feature. This means that, while integrating the Initial Value Problem
:eq:`CVODE_ivp`, CVODE can also find the roots of a set of
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
CVODE. If such a root is desired, the user should reformulate the
root function so that it changes sign at the desired root.

The basic scheme used is to check for sign changes of any :math:`g_i(t)`
over each time step taken, and then (when a sign change is found) to
hone in on the root(s) with a modified secant method
:cite:p:`HeSh:80`. In addition, each time :math:`g` is
computed, CVODE checks to see if :math:`g_i(t) = 0` exactly, and if
so it reports this as a root. However, if an exact zero of any
:math:`g_i` is found at a point :math:`t`, CVODE computes :math:`g`
at :math:`t + \delta` for a small increment :math:`\delta`, slightly
further in the direction of integration, and if any
:math:`g_i(t + \delta)=0` also, CVODE stops and reports an error.
This way, each time CVODE takes a time step, it is guaranteed that
the values of all :math:`g_i` are nonzero at some past value of
:math:`t`, beyond which a search for roots is to be done.

At any given time in the course of the time-stepping, after suitable
checking and adjusting has been done, CVODE has an interval
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
