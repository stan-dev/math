.. ----------------------------------------------------------------
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

.. _ARKODE.Mathematics:

===========================
Mathematical Considerations
===========================

ARKODE solves ODE initial value problems (IVP) in :math:`\mathbb{R}^N`
posed in the form

.. math::
   M(t)\, \dot{y} = f(t,y), \qquad y(t_0) = y_0.
   :label: ARKODE_IVP

Here, :math:`t` is the independent variable (e.g. time), and the
dependent variables are given by :math:`y \in \mathbb{R}^N`, where we
use the notation :math:`\dot{y}` to denote :math:`\mathrm dy/\mathrm dt`.

For each value of :math:`t`, :math:`M(t)` is a user-specified linear
operator from :math:`\mathbb{R}^N \to \mathbb{R}^N`.  This operator
is assumed to be nonsingular and independent of :math:`y`.  For
standard systems of ordinary differential equations and for
problems arising from the spatial semi-discretization of partial
differential equations using finite difference, finite volume, or
spectral finite element methods, :math:`M` is typically the identity
matrix, :math:`I`.  For PDEs using standard finite-element
spatial semi-discretizations, :math:`M` is typically a
well-conditioned mass matrix that is fixed throughout a simulation
(or at least fixed between spatial rediscretization events).

The ODE right-hand side is given by the function :math:`f(t,y)` --
in general we make no assumption that the problem :eq:`ARKODE_IVP` is
autonomous (i.e., :math:`f=f(y)`) or linear (:math:`f=Ay`).
In general, the time integration methods within ARKODE support
additive splittings of this right-hand side function, as described
in the subsections that follow.  Through these splittings, the
time-stepping methods currently supplied with ARKODE are designed
to solve stiff, nonstiff, mixed stiff/nonstiff, and multirate
problems.  As per Ascher and Petzold :cite:p:`AsPe:98`, a problem is "stiff"
if the stepsize needed to maintain stability of the forward Euler
method is much smaller than that required to represent the solution
accurately.

In the sub-sections that follow, we elaborate on the numerical
methods utilized in ARKODE.  We first discuss the "single-step" nature
of the ARKODE infrastructure, including its usage modes and approaches
for interpolated solution output.  We then discuss the current suite
of time-stepping modules supplied with ARKODE, including the ARKStep
module for :ref:`additive Runge--Kutta methods <ARKODE.Mathematics.ARK>`,
the ERKStep module that is optimized for :ref:`explicit Runge--Kutta
methods <ARKODE.Mathematics.ERK>`, and the MRIStep module for :ref:`multirate
infinitesimal step (MIS), multirate infinitesimal GARK (MRI-GARK), and
implicit-explicit MRI-GARK (IMEX-MRI-GARK) methods <ARKODE.Mathematics.MRIStep>`.
We then discuss the :ref:`adaptive temporal error controllers
<ARKODE.Mathematics.Adaptivity>` shared by the time-stepping modules, including
discussion of our choice of norms for measuring errors within various components
of the solver.

We then discuss the nonlinear and linear solver strategies used by
ARKODE's time-stepping modules for solving implicit algebraic systems
that arise in computing each stage and/or step:
:ref:`nonlinear solvers <ARKODE.Mathematics.Nonlinear>`,
:ref:`linear solvers <ARKODE.Mathematics.Linear>`,
:ref:`preconditioners <ARKODE.Mathematics.Preconditioning>`,
:ref:`error control <ARKODE.Mathematics.Error>` within iterative nonlinear
and linear solvers, algorithms for
:ref:`initial predictors <ARKODE.Mathematics.Predictors>` for implicit stage
solutions, and approaches for handling
:ref:`non-identity mass-matrices <ARKODE.Mathematics.MassSolve>`.

We conclude with a section describing ARKODE's :ref:`rootfinding
capabilities <ARKODE.Mathematics.Rootfinding>`, that may be used to stop
integration of a problem prematurely based on traversal of roots in
user-specified functions.



.. _ARKODE.Mathematics.SingleStep:

Adaptive single-step methods
===============================

The ARKODE infrastructure is designed to support single-step, IVP
integration methods, i.e.

.. math::

   y_{n} = \varphi(y_{n-1}, h_n)

where :math:`y_{n-1}` is an approximation to the solution :math:`y(t_{n-1})`,
:math:`y_{n}` is an approximation to the solution :math:`y(t_n)`,
:math:`t_n = t_{n-1} + h_n`, and the approximation method is
represented by the function :math:`\varphi`.

The choice of step size :math:`h_n` is determined by the time-stepping
method (based on user-provided inputs, typically accuracy requirements).
However, users may place minimum/maximum bounds on :math:`h_n` if desired.

ARKODE's time stepping modules may be run in a variety of "modes":

* **NORMAL** -- The solver will take internal steps until it has just
  overtaken a user-specified output time, :math:`t_\text{out}`, in the
  direction of integration, i.e. :math:`t_{n-1} < t_\text{out} \le
  t_{n}` for forward integration, or :math:`t_{n} \le t_\text{out} <
  t_{n-1}` for backward integration.  It will then compute an
  approximation to the solution :math:`y(t_\text{out})` by
  interpolation (using one of the dense output routines described in
  the section :numref:`ARKODE.Mathematics.Interpolation`).

* **ONE-STEP** -- The solver will only take a single internal step
  :math:`y_{n-1} \to y_{n}` and then return control back to the
  calling program.  If this step will overtake :math:`t_\text{out}`
  then the solver will again return an interpolated result; otherwise
  it will return a copy of the internal solution :math:`y_{n}`.

* **NORMAL-TSTOP** -- The solver will take internal steps until the next
  step will overtake :math:`t_\text{out}`.  It will then limit
  this next step so that :math:`t_n = t_{n-1} + h_n = t_\text{out}`,
  and once the step completes it will return a copy of the internal
  solution :math:`y_{n}`.

* **ONE-STEP-TSTOP** -- The solver will check whether the next step
  will overtake :math:`t_\text{out}` -- if not then this mode is
  identical to "one-step" above; otherwise it will limit this next
  step so that :math:`t_n = t_{n-1} + h_n = t_\text{out}`.  In either
  case, once the step completes it will return a copy of the internal
  solution :math:`y_{n}`.

We note that interpolated solutions may be slightly less accurate than
the internal solutions produced by the solver.  Hence, to ensure that
the returned value has full method accuracy one of the "tstop" modes
may be used.



.. _ARKODE.Mathematics.Interpolation:

Interpolation
===============

As mentioned above, the time-stepping modules in ARKODE support
interpolation of solutions :math:`y(t_\text{out})` and derivatives
:math:`y^{(d)}(t_\text{out})`, where :math:`t_\text{out}` occurs
within a completed time step from :math:`t_{n-1} \to t_n`.
Additionally, this module supports extrapolation of solutions and
derivatives for :math:`t` outside this interval (e.g. to construct
predictors for iterative nonlinear and linear solvers).  To this end,
ARKODE currently supports construction of polynomial interpolants
:math:`p_q(t)` of polynomial degree up to :math:`q=5`, although
users may select interpolants of lower degree.

ARKODE provides two complementary interpolation approaches,
both of which are accessible from any of the
time-stepping modules: "Hermite" and "Lagrange".  The former approach
has been included with ARKODE since its inception, and is more
suitable for non-stiff problems; the latter is a new approach that is
designed to provide increased accuracy when integrating stiff problems.
Both are described in detail below.


.. _ARKODE.Mathematics.Interpolation.Hermite:

Hermite interpolation module
-----------------------------

For non-stiff problems, polynomial interpolants of Hermite form are provided.
Rewriting the IVP :eq:`ARKODE_IVP` in standard form,

.. math::
   \dot{y} = \hat{f}(t,y), \qquad y(t_0) = y_0.

we typically construct temporal interpolants using the data
:math:`\left\{ y_{n-1}, \hat{f}_{n-1}, y_{n}, \hat{f}_{n} \right\}`,
where here we use the simplified notation :math:`\hat{f}_{k}` to denote
:math:`\hat{f}(t_k,y_k)`.  Defining a normalized "time" variable,
:math:`\tau`, for the most-recently-computed solution interval
:math:`t_{n-1} \to t_{n}` as

.. math::

   \tau(t) = \frac{t-t_{n}}{h_{n}},

we then construct the interpolants :math:`p_q(t)` as follows:

* :math:`q=0`: constant interpolant

  .. math::

     p_0(\tau) = \frac{y_{n-1} + y_{n}}{2}.

* :math:`q=1`: linear Lagrange interpolant

  .. math::

     p_1(\tau) = -\tau\, y_{n-1} + (1+\tau)\, y_{n}.

* :math:`q=2`: quadratic Hermite interpolant

  .. math::

     p_2(\tau) =  \tau^2\,y_{n-1} + (1-\tau^2)\,y_{n} + h_n(\tau+\tau^2)\,\hat{f}_{n}.

* :math:`q=3`: cubic Hermite interpolant

  .. math::

     p_3(\tau) =  (3\tau^2 + 2\tau^3)\,y_{n-1} +
     (1-3\tau^2-2\tau^3)\,y_{n} + h_n(\tau^2+\tau^3)\,\hat{f}_{n-1} +
     h_n(\tau+2\tau^2+\tau^3)\,\hat{f}_{n}.

* :math:`q=4`: quartic Hermite interpolant

  .. math::

     p_4(\tau) &= (-6\tau^2 - 16\tau^3 - 9\tau^4)\,y_{n-1} +
     (1 + 6\tau^2 + 16\tau^3 + 9\tau^4)\,y_{n} +
     \frac{h_n}{4}(-5\tau^2 - 14\tau^3 - 9\tau^4)\,\hat{f}_{n-1} \\
     &+ h_n(\tau + 2\tau^2 + \tau^3)\,\hat{f}_{n} +
     \frac{27 h_n}{4}(-\tau^4 - 2\tau^3 - \tau^2)\,\hat{f}_a,

  where :math:`\hat{f}_a=\hat{f}\left(t_{n} - \dfrac{h_n}{3},p_3\left(-\dfrac13\right)\right)`.
  We point out that interpolation at this degree requires an additional evaluation
  of the full right-hand side function :math:`\hat{f}(t,y)`, thereby increasing its
  cost in comparison with :math:`p_3(t)`.

* :math:`q=5`: quintic Hermite interpolant

  .. math::

     p_5(\tau) &= (54\tau^5 + 135\tau^4 + 110\tau^3 + 30\tau^2)\,y_{n-1} +
     (1 - 54\tau^5 - 135\tau^4 - 110\tau^3 - 30\tau^2)\,y_{n} \\
     &+ \frac{h_n}{4}(27\tau^5 + 63\tau^4 + 49\tau^3 + 13\tau^2)\,\hat{f}_{n-1} +
     \frac{h_n}{4}(27\tau^5 + 72\tau^4 + 67\tau^3 + 26\tau^2 + \tau)\,\hat{f}_n \\
     &+ \frac{h_n}{4}(81\tau^5 + 189\tau^4 + 135\tau^3 + 27\tau^2)\,\hat{f}_a +
     \frac{h_n}{4}(81\tau^5 + 216\tau^4 + 189\tau^3 + 54\tau^2)\,\hat{f}_b,

  where :math:`\hat{f}_a=\hat{f}\left(t_{n} - \dfrac{h_n}{3},p_4\left(-\dfrac13\right)\right)`
  and :math:`\hat{f}_b=\hat{f}\left(t_{n} - \dfrac{2h_n}{3},p_4\left(-\dfrac23\right)\right)`.
  We point out that interpolation at this degree requires four additional evaluations
  of the full right-hand side function :math:`\hat{f}(t,y)`, thereby significantly
  increasing its cost over :math:`p_4(t)`.

We note that although interpolants of order :math:`q > 5` are possible, these are
not currently implemented due to their increased computing and storage costs.



.. _ARKODE.Mathematics.Interpolation.Lagrange:

Lagrange interpolation module
-----------------------------

For stiff problems where :math:`\hat{f}` may have large Lipschitz constant,
polynomial interpolants of Lagrange form are provided.  These interpolants
are constructed using the data
:math:`\left\{ y_{n}, y_{n-1}, \ldots, y_{n-\nu} \right\}` where
:math:`0\le\nu\le5`.  These polynomials have the form

.. math::

   p(t) &= \sum_{j=0}^{\nu} y_{n-j} p_j(t),\quad\text{where}\\
   p_j(t) &= \prod_{\substack{l=0\\ l\ne j}}^{\nu} \left(\frac{t-t_l}{t_j-t_l}\right), \quad j=0,\ldots,\nu.

Since we assume that the solutions :math:`y_{n-j}` have length much larger
than :math:`\nu\le5` in ARKODE-based simulations, we evaluate :math:`p` at
any desired :math:`t\in\mathbb{R}` by first evaluating the Lagrange polynomial
basis functions at the input value for :math:`t`, and then performing a simple linear
combination of the vectors :math:`\{y_k\}_{k=0}^{\nu}`.  Derivatives :math:`p^{(d)}(t)`
may be evaluated similarly as

.. math::

   p^{(d)}(t) = \sum_{j=0}^{\nu} y_{n-j}\, p_j^{(d)}(t),

however since the algorithmic complexity involved in evaluating derivatives of the
Lagrange basis functions increases dramatically as the derivative order grows, our Lagrange
interpolation module currently only provides derivatives up to :math:`d=3`.

We note that when using this interpolation module, during the first
:math:`(\nu-1)` steps of integration we do not have sufficient solution history
to construct the full :math:`\nu`-degree interpolant.  Therefore during these
initial steps, we construct the highest-degree interpolants that are currently
available at the moment, achieving the full :math:`\nu`-degree interpolant once
these initial steps have completed.



.. _ARKODE.Mathematics.ARK:

ARKStep -- Additive Runge--Kutta methods
=========================================

The ARKStep time-stepping module in ARKODE is designed for IVPs of the
form

.. math::
   M(t)\, \dot{y} = f^E(t,y) + f^I(t,y), \qquad y(t_0) = y_0,
   :label: ARKODE_IMEX_IVP

i.e. the right-hand side function is additively split into two
components:

* :math:`f^E(t,y)` contains the "nonstiff" components of the
  system (this will be integrated using an explicit method);

* :math:`f^I(t,y)` contains the "stiff" components of the
  system (this will be integrated using an implicit method);

and the left-hand side may include a nonsingular, possibly
time-dependent,  matrix :math:`M(t)`.

In solving the IVP :eq:`ARKODE_IMEX_IVP`, we first consider the corresponding
problem in standard form,

.. math::
   \dot{y} = \hat{f}^E(t,y) + \hat{f}^I(t,y), \qquad y(t_0) = y_0,
   :label: ARKODE_IMEX_IVP_standard

where :math:`\hat{f}^E(t,y) = M(t)^{-1}\,f^E(t,y)` and
:math:`\hat{f}^I(t,y) = M(t)^{-1}\,f^I(t,y)`.  ARKStep then utilizes variable-step,
embedded, :index:`additive Runge--Kutta methods` (ARK), corresponding
to algorithms of the form

.. math::
   z_i &= y_{n-1} + h_n \sum_{j=1}^{i-1} A^E_{i,j} \hat{f}^E(t^E_{n,j}, z_j)
                  + h_n \sum_{j=1}^{i} A^I_{i,j} \hat{f}^I(t^I_{n,j}, z_j),
   \quad i=1,\ldots,s, \\
   y_n &= y_{n-1} + h_n \sum_{i=1}^{s} \left(b^E_i \hat{f}^E(t^E_{n,i}, z_i)
                 + b^I_i \hat{f}^I(t^I_{n,i}, z_i)\right), \\
   \tilde{y}_n &= y_{n-1} + h_n \sum_{i=1}^{s} \left(
                  \tilde{b}^E_i \hat{f}^E(t^E_{n,i}, z_i) +
                  \tilde{b}^I_i \hat{f}^I(t^I_{n,i}, z_i)\right).
   :label: ARKODE_ARK

Here :math:`\tilde{y}_n` are embedded solutions that approximate
:math:`y(t_n)` and are used for error estimation; these typically
have slightly lower accuracy than the computed solutions :math:`y_n`.
The internal stage times are abbreviated using the notation
:math:`t^E_{n,j} = t_{n-1} + c^E_j h_n` and
:math:`t^I_{n,j} = t_{n-1} + c^I_j h_n`.  The ARK method is
primarily defined through the coefficients :math:`A^E \in
\mathbb{R}^{s\times s}`, :math:`A^I \in \mathbb{R}^{s\times s}`,
:math:`b^E \in \mathbb{R}^{s}`, :math:`b^I \in \mathbb{R}^{s}`,
:math:`c^E \in \mathbb{R}^{s}` and :math:`c^I \in \mathbb{R}^{s}`,
that correspond with the explicit and implicit Butcher tables.
Additional coefficients :math:`\tilde{b}^E \in \mathbb{R}^{s}` and
:math:`\tilde{b}^I \in \mathbb{R}^{s}` are used to construct the
embedding :math:`\tilde{y}_n`.  We note that ARKStep currently
enforces the constraint that the explicit and implicit methods in an
ARK pair must share the same number of stages, :math:`s`.  We note that
except when the problem has a time-independent mass matrix :math:`M`, ARKStep
allows the possibility for different explicit and implicit abscissae,
i.e. :math:`c^E` need not equal :math:`c^I`.

The user of ARKStep must choose appropriately between one of three
classes of methods: *ImEx*, *explicit*, and *implicit*.  All of
the built-in Butcher tables encoding the coefficients
:math:`c^E`, :math:`c^I`, :math:`A^E`, :math:`A^I`, :math:`b^E`,
:math:`b^I`, :math:`\tilde{b}^E` and :math:`\tilde{b}^I` are further
described in the section :numref:`Butcher`.

For mixed stiff/nonstiff problems, a user should provide both of the
functions :math:`f^E` and :math:`f^I` that define the IVP system.  For
such problems, ARKStep currently implements the ARK methods proposed in
:cite:p:`KenCarp:03`, allowing for methods having order of accuracy :math:`q =
\{3,4,5\}` and embeddings with orders :math:`p = \{2, 3, 4\}`;
the tables for these methods are given in section :numref:`Butcher.additive`.
Additionally, user-defined ARK tables are supported.

For nonstiff problems, a user may specify that :math:`f^I = 0`,
i.e. the equation :eq:`ARKODE_IMEX_IVP` reduces to the non-split IVP

.. math::
   M(t)\, \dot{y} = f^E(t,y), \qquad y(t_0) = y_0.
   :label: ARKODE_IVP_explicit

In this scenario, the coefficients :math:`A^I=0`, :math:`c^I=0`,
:math:`b^I=0` and :math:`\tilde{b}^I=0` in :eq:`ARKODE_ARK`, and the ARK
methods reduce to classical :index:`explicit Runge--Kutta methods`
(ERK).  For these classes of methods, ARKODE provides coefficients
with orders of accuracy :math:`q = \{2,3,4,5,6,8\}`, with embeddings
of orders :math:`p = \{1,2,3,4,5,7\}`.  These default to the methods in
sections
:numref:`Butcher.Heun_Euler`,
:numref:`Butcher.Bogacki_Shampine`, :numref:`Butcher.Zonneveld`,
:numref:`Butcher.Cash-Karp`, :numref:`Butcher.Verner-6-5`, and
:numref:`Butcher.Fehlberg-8-7`, respectively.  As with ARK
methods, user-defined ERK tables are supported.

Alternately, for stiff problems the user may specify that :math:`f^E = 0`,
so the equation :eq:`ARKODE_IMEX_IVP` reduces to the non-split IVP

.. math::
   M(t)\, \dot{y} = f^I(t,y), \qquad y(t_0) = y_0.
   :label: ARKODE_IVP_implicit

Similarly to ERK methods, in this scenario the coefficients
:math:`A^E=0`, :math:`c^E=0`, :math:`b^E=0` and :math:`\tilde{b}^E=0`
in :eq:`ARKODE_ARK`, and the ARK methods reduce to classical
:index:`diagonally-implicit Runge--Kutta methods` (DIRK).  For these
classes of methods, ARKODE provides tables with orders of accuracy
:math:`q = \{2,3,4,5\}`, with embeddings of orders
:math:`p = \{1,2,3,4\}`. These default to the methods
:numref:`Butcher.SDIRK-2-1`, :numref:`Butcher.ARK_4_2_3_I`,
:numref:`Butcher.SDIRK-5-4`, and :numref:`Butcher.ARK_8_4_5_I`,
respectively.  Again, user-defined DIRK tables are supported.



.. _ARKODE.Mathematics.ERK:

ERKStep -- Explicit Runge--Kutta methods
===========================================

The ERKStep time-stepping module in ARKODE is designed for IVP
of the form

.. math::
   \dot{y} = f(t,y), \qquad y(t_0) = y_0,
   :label: ARKODE_IVP_simple_explicit

i.e., unlike the more general problem form :eq:`ARKODE_IMEX_IVP`, ERKStep
requires that problems have an identity mass matrix (i.e., :math:`M(t)=I`)
and that the right-hand side function is not split into separate
components.

For such problems, ERKStep provides variable-step, embedded,
:index:`explicit Runge--Kutta methods` (ERK), corresponding to
algorithms of the form

.. math::
   z_i &= y_{n-1} + h_n \sum_{j=1}^{i-1} A_{i,j} f(t_{n,j}, z_j),
   \quad i=1,\ldots,s, \\
   y_n &= y_{n-1} + h_n \sum_{i=1}^{s} b_i f(t_{n,i}, z_i), \\
   \tilde{y}_n &= y_{n-1} + h_n \sum_{i=1}^{s} \tilde{b}_i f(t_{n,i}, z_i),
   :label: ARKODE_ERK

where the variables have the same meanings as in the previous section.

Clearly, the problem :eq:`ARKODE_IVP_simple_explicit` is fully encapsulated
in the more general problem :eq:`ARKODE_IVP_explicit`, and the algorithm
:eq:`ARKODE_ERK` is similarly encapsulated in the more general algorithm :eq:`ARKODE_ARK`.
While it therefore follows that ARKStep can be used to solve every
problem solvable by ERKStep, using the same set of methods, we
include ERKStep as a distinct time-stepping module since this
simplified form admits a more efficient and memory-friendly implementation
than the more general form :eq:`ARKODE_IVP_simple_explicit`.


.. _ARKODE.Mathematics.MRIStep:

MRIStep -- Multirate infinitesimal step methods
================================================

The MRIStep time-stepping module in ARKODE is designed for IVPs
of the form

.. math::
   \dot{y} = f^E(t,y) + f^I(t,y) + f^F(t,y), \qquad y(t_0) = y_0.
   :label: ARKODE_IVP_two_rate

i.e., the right-hand side function is additively split into three
components:

* :math:`f^E(t,y)` contains the "slow-nonstiff" components of the system
  (this will be integrated using an explicit method and a large time step
  :math:`h^S`),

* :math:`f^I(t,y)` contains the "slow-stiff" components of the system
  (this will be integrated using an implicit method and a large time step
  :math:`h^S`), and

* :math:`f^F(t,y)` contains the "fast" components of the system (this will be
  integrated using a possibly different method than the slow time scale and a
  small time step :math:`h^F \ll h^S`).

As with ERKStep, MRIStep currently requires that problems be posed with
an identity mass matrix, :math:`M(t)=I`. The slow time scale may consist of only
nonstiff terms (:math:`f^I \equiv 0`), only stiff terms (:math:`f^E \equiv 0`),
or both nonstiff and stiff terms.

For cases with only a single slow right-hand side function (i.e.,
:math:`f^E \equiv 0` or :math:`f^I \equiv 0`), MRIStep provides fixed-slow-step
multirate infinitesimal step (MIS) :cite:p:`Schlegel:09, Schlegel:12a,
Schlegel:12b` and multirate infinitesimal GARK (MRI-GARK) :cite:p:`Sandu:19`
methods. For problems with an additively split slow right-hand side MRIStep
provides fixed-slow-step implicit-explicit MRI-GARK (IMEX-MRI-GARK)
:cite:p:`ChiRen:21` methods.  The slow (outer) method derives from an :math:`s`
stage Runge--Kutta method for MIS and MRI-GARK methods or an additive Runge--Kutta
method for IMEX-MRI-GARK methods. In either case, the stage values and the new
solution are computed by solving an auxiliary ODE with a fast (inner) time
integration method. This corresponds to the following algorithm for a single
step:

#. Set :math:`z_1 = y_{n-1}`.

#. For :math:`i = 2,\ldots,s+1` do:

   #. Let :math:`t_{n,i-1}^S = t_{n-1} + c_{i-1}^S h^S` and
      :math:`v(t_{n,i-1}^S) = z_{i-1}`.

   #. Let :math:`r_i(t) =
      \frac{1}{\Delta c_i^S}
      \sum\limits_{j=1}^{i-1} \omega_{i,j}(\tau) f^E(t_{n,j}^I, z_j) +
      \frac{1}{\Delta c_i^S}
      \sum\limits_{j=1}^i \gamma_{i,j}(\tau) f^I(t_{n,j}^I, z_j)`
      where :math:`\Delta c_i^S=\left(c^S_i - c^S_{i-1}\right)` and the
      normalized time is :math:`\tau = (t - t_{n,i-1}^S)/(h^S \Delta c_i^S)`.

   #. For :math:`t \in [t_{n,i-1}^S, t_{n,i}^S]` solve
      :math:`\dot{v}(t) = f^F(t, v) + r_i(t)`.

   #. Set :math:`z_i = v(t_{n,i}^S)`.

#. Set :math:`y_{n} = z_{s+1}`.

The fast (inner) IVP solve can be carried out using either the ARKStep module
(allowing for explicit, implicit, or ImEx treatments of the fast time scale with
fixed or adaptive steps), or a user-defined integration method (see section
:numref:`ARKODE.Usage.MRIStep.CustomInnerStepper`).

The final abscissa is :math:`c^S_{s+1}=1` and the coefficients
:math:`\omega_{i,j}` and :math:`\gamma_{i,j}` are polynomials in time that
dictate the couplings from the slow to the fast time scale; these can be
expressed as in :cite:p:`ChiRen:21` and :cite:p:`Sandu:19` as

.. math::
   \omega_{i,j}(\tau) = \sum_{k\geq 0} \omega_{i,j}^{\{k\}} \tau^k
   \quad\text{and}\quad
   \gamma_{i,j}(\tau) = \sum_{k\geq 0} \gamma_{i,j}^{\{k\}} \tau^k,
   :label: ARKODE_MRI_coupling

and where the tables :math:`\Omega^{\{k\}}\in\mathbb{R}^{(s+1)\times(s+1)}` and
:math:`\Gamma^{\{k\}}\in\mathbb{R}^{(s+1)\times(s+1)}` define the slow-to-fast
coupling for the explicit and implicit components respectively.

For traditional MIS methods, the coupling coefficients are uniquely defined
based on a slow Butcher table :math:`(A^S,b^S,c^S)` having an explicit first
stage (i.e., :math:`c^S_1=0` and :math:`A^S_{1,j}=0` for :math:`1\le j\le s`),
sorted abscissae (i.e., :math:`c^S_{i} \ge c^S_{i-1}` for :math:`2\le i\le s`),
and the final abscissa is :math:`c^S_s \leq 1`. With these properties met, the
coupling coefficients for an explicit-slow method are given as

.. math::
   \omega_{i,j}^{\{0\}} = \begin{cases}
   0, & \text{if}\; i=1,\\
   A^S_{i,j} - A^S_{i-1,j}, & \text{if}\; 2\le i\le s,\\
   b^S_j - A^S_{s,j}, & \text{if}\; i=s+1.
   \end{cases}
   :label: ARKODE_MIS_to_MRI

For general slow tables :math:`(A^S,b^S,c^S)` with at least second-order
accuracy, the corresponding MIS method will be second order. However, if this
slow table is at least third order and satisfies the additional condition

.. math::
   \sum_{i=2}^{s} \left(c_i^S-c_{i-1}^S\right)
   \left(\mathbf{e}_i+\mathbf{e}_{i-1}\right)^T A^S c^S
   + \left(1-c_{s}^S\right) \left( \frac12+\mathbf{e}_{s}^T A^S c^S \right)
   = \frac13,
   :label: ARKODE_MIS_order3

where :math:`\mathbf{e}_j` corresponds to the :math:`j`-th column from the
:math:`s \times s` identity matrix, then the overall MIS method will be third
order.

In the above algorithm, when the slow (outer) method has repeated abscissa, i.e.
:math:`\Delta c_i^S = 0` for stage :math:`i`, the fast (inner) IVP can be
rescaled and integrated analytically. In this case the stage is computed as

.. math::
   z_i = z_{i-1}
   + h^S \sum_{j=1}^{i-1} \left(\sum_{k\geq 0}
     \frac{\omega_{i,j}^{\{k\}}}{k+1}\right) f^E(t_{n,j}^S, z_j)
   + h^S \sum_{j=1}^i \left(\sum_{k\geq 0}
     \frac{\gamma_{i,j}^{\{k\}}}{k+1}\right) f^I(t_{n,j}^S, z_j),
   :label: ARKODE_MRI_delta_c_zero

which corresponds to a standard ARK, DIRK, or ERK stage computation depending on
whether the summations over :math:`k` are zero or nonzero.

As with standard ARK and DIRK methods, implicitness at the slow time scale is
characterized by nonzero values on or above the diagonal of the matrices
:math:`\Gamma^{\{k\}}`. Typically, MRI-GARK and IMEX-MRI-GARK methods are at
most diagonally-implicit (i.e., :math:`\gamma_{i,j}^{\{k\}}=0` for all
:math:`j>i`). Furthermore, diagonally-implicit stages are characterized as being
"solve-decoupled" if :math:`\Delta c_i^S = 0` when `\gamma_{i,i}^{\{k\}} \ne 0`,
in which case the stage is computed as standard ARK or DIRK update. Alternately,
a diagonally-implicit stage :math:`i` is considered "solve-coupled" if
:math:`\Delta c^S_i \gamma_{i,j}^{\{k\}} \ne 0`, in which
case the stage solution :math:`z_i` is *both* an input to :math:`r(t)` and the
result of time-evolution of the fast IVP, necessitating an implicit solve that
is coupled to the fast (inner) solver. At present, only "solve-decoupled"
diagonally-implicit MRI-GARK and IMEX-MRI-GARK methods are supported.

For problems with only a slow-nonstiff term (:math:`f^I \equiv 0`), MRIStep
provides third and fourth order explicit MRI-GARK methods. In cases with only a
slow-stiff term (:math:`f^E \equiv 0`), MRIStep supplies second, third, and
fourth order implicit solve-decoupled MRI-GARK methods. For applications
with both stiff and nonstiff slow terms, MRIStep implements third and fourth
order IMEX-MRI-GARK methods. For a complete list of the methods available in
MRIStep see :numref:`ARKODE.Usage.MRIStep.MRIStepCoupling.Tables`. Additionally, users
may supply their own method by defining and attaching a coupling table, see
:numref:`ARKODE.Usage.MRIStep.MRIStepCoupling` for more information.


.. _ARKODE.Mathematics.Error.Norm:

Error norms
============================

In the process of controlling errors at various levels (time
integration, nonlinear solution, linear solution), the methods in
ARKODE use a :index:`weighted root-mean-square norm`, denoted
:math:`\|\cdot\|_\text{WRMS}`, for all error-like quantities,

.. math::
   \|v\|_\text{WRMS} = \left( \frac{1}{N} \sum_{i=1}^N \left(v_i\,
   w_i\right)^2\right)^{1/2}.
   :label: ARKODE_WRMS_NORM

The utility of this norm arises in the specification of the weighting
vector :math:`w`, that combines the units of the problem with
user-supplied values that specify an "acceptable" level of error.  To
this end, we construct an :index:`error weight vector` using
the most-recent step solution and user-supplied relative and
absolute tolerances, namely

.. math::
   w_i = \big(RTOL\cdot |y_{n-1,i}| + ATOL_i\big)^{-1}.
   :label: ARKODE_EWT

Since :math:`1/w_i` represents a tolerance in the :math:`i`-th component of the
solution vector :math:`y`, a vector whose WRMS norm is 1 is regarded
as "small."  For brevity, unless specified otherwise we will drop the
subscript WRMS on norms in the remainder of this section.

Additionally, for problems involving a non-identity mass matrix,
:math:`M\ne I`, the units of equation :eq:`ARKODE_IMEX_IVP` may differ from the
units of the solution :math:`y`.  In this case, we may additionally
construct a :index:`residual weight vector`,

.. math::
   w_i = \Big(RTOL\cdot \left| \big(M(t_{n-1})\, y_{n-1}\big)_i \right| + ATOL'_i\Big)^{-1},
   :label: ARKODE_RWT

where the user may specify a separate absolute residual tolerance
value or array, :math:`ATOL'`.  The choice of weighting vector used
in any given norm is determined by the quantity being measured: values
having "solution" units use :eq:`ARKODE_EWT`, whereas values having "equation"
units use :eq:`ARKODE_RWT`.  Obviously, for problems with :math:`M=I`, the
solution and equation units are identical, in which case the solvers in
ARKODE will use :eq:`ARKODE_EWT` when computing all error norms.




.. _ARKODE.Mathematics.Adaptivity:

Time step adaptivity
=======================

A critical component of IVP "solvers" (rather than just
time-steppers) is their adaptive control of local truncation error (LTE).
At every step, we estimate the local error, and ensure that it
satisfies tolerance conditions.  If this local error test fails, then
the step is recomputed with a reduced step size.  To this end, the
Runge--Kutta methods packaged within both the ARKStep and ERKStep
modules admit an embedded solution :math:`\tilde{y}_n`, as shown in
equations :eq:`ARKODE_ARK` and :eq:`ARKODE_ERK`.  Generally, these embedded
solutions attain a slightly lower order of accuracy than the computed
solution :math:`y_n`.  Denoting the order of accuracy for :math:`y_n`
as :math:`q` and for :math:`\tilde{y}_n` as :math:`p`, most of these
embedded methods satisfy :math:`p = q-1`.  These values of :math:`q`
and :math:`p` correspond to the *global* orders of accuracy for the
method and embedding, hence each admit local truncation errors
satisfying :cite:p:`HWN:87`

.. math::
   \| y_n - y(t_n) \| = C h_n^{q+1} + \mathcal O(h_n^{q+2}), \\
   \| \tilde{y}_n - y(t_n) \| = D h_n^{p+1} + \mathcal O(h_n^{p+2}),
   :label: ARKODE_AsymptoticErrors

where :math:`C` and :math:`D` are constants independent of
:math:`h_n`, and where we have assumed exact initial conditions for
the step, i.e. :math:`y_{n-1} = y(t_{n-1})`. Combining these
estimates, we have

.. math::
   \| y_n - \tilde{y}_n \| = \| y_n - y(t_n) - \tilde{y}_n + y(t_n) \|
   \le \| y_n - y(t_n) \| + \| \tilde{y}_n - y(t_n) \|
   \le D h_n^{p+1} + \mathcal O(h_n^{p+2}).

We therefore use the norm of the difference between :math:`y_n` and
:math:`\tilde{y}_n` as an estimate for the LTE at the step :math:`n`

.. math::
   T_n = \beta \left(y_n - \tilde{y}_n\right) =
   \beta h_n \sum_{i=1}^{s} \left[
   \left(b^E_i - \tilde{b}^E_i\right) \hat{f}^E(t^E_{n,i}, z_i) +
   \left(b^I_i - \tilde{b}^I_i\right) \hat{f}^I(t^I_{n,i}, z_i) \right]
   :label: ARKODE_LTE

for ARK methods, and similarly for ERK methods.  Here, :math:`\beta>0`
is an error *bias* to help account for the error constant :math:`D`;
the default value of this constant is :math:`\beta = 1.5`, which may
be modified by the user.

With this LTE estimate, the local error test is simply
:math:`\|T_n\| < 1` since this norm includes the user-specified
tolerances.  If this error test passes, the step is considered
successful, and the estimate is subsequently used to determine the next
step size, the algorithms used for this purpose are described in
:numref:`ARKODE.Mathematics.Adaptivity`.  If the error
test fails, the step is rejected and a new step size :math:`h'` is
then computed using the same error controller as for successful steps.
A new attempt at the step is made, and the error test is repeated.  If
the error test fails twice, then :math:`h'/h` is limited above to 0.3,
and limited below to 0.1 after an additional step failure.  After
seven error test failures, control is returned to the user with a
failure message.  We note that all of the constants listed above are
only the default values; each may be modified by the user.

We define the step size ratio between a prospective step :math:`h'`
and a completed step :math:`h` as :math:`\eta`, i.e. :math:`\eta = h'
/ h`.  This value is subsequently bounded from above by
:math:`\eta_\text{max}` to ensure that step size adjustments are not
overly aggressive.  This upper bound changes according to the step
and history,

.. math::
   \eta_\text{max} = \begin{cases}
     \text{etamx1}, & \quad\text{on the first step (default is 10000)}, \\
     \text{growth}, & \quad\text{on general steps (default is 20)}, \\
     1, & \quad\text{if the previous step had an error test failure}.
   \end{cases}

A flowchart detailing how the time steps are modified at each
iteration to ensure solver convergence and successful steps is given
in the figure below.  Here, all norms correspond to the WRMS norm, and
the error adaptivity function **arkAdapt** is supplied by one of the
error control algorithms discussed in the subsections below.

.. _adaptivity_figure:
.. figure:: /figs/arkode/time_adaptivity.png
   :width: 50%
   :align: center


For some problems it may be preferable to avoid small step size
adjustments.  This can be especially true for problems that construct
a Newton Jacobian matrix or a preconditioner for a nonlinear or an
iterative linear solve, where this construction is computationally
expensive, and where convergence can be seriously hindered through use
of an inaccurate matrix.  To accommodate these scenarios, the step is
left unchanged when :math:`\eta \in [\eta_L, \eta_U]`.  The default
values for this interval are :math:`\eta_L = 1` and :math:`\eta_U =
1.5`, and may be modified by the user.

We note that any choices for :math:`\eta` (or equivalently,
:math:`h'`) are subsequently constrained by the optional user-supplied
bounds :math:`h_\text{min}` and :math:`h_\text{max}`.  Additionally,
the time-stepping algorithms in ARKODE may similarly limit :math:`h'`
to adhere to a user-provided "TSTOP" stopping point,
:math:`t_\text{stop}`.



The time-stepping modules in ARKODE adapt the step
size in order to attain local errors within desired tolerances of the
true solution.  These adaptivity algorithms estimate the prospective
step size :math:`h'` based on the asymptotic local error estimates
:eq:`ARKODE_AsymptoticErrors`.  We define the values :math:`\varepsilon_n`,
:math:`\varepsilon_{n-1}` and :math:`\varepsilon_{n-2}` as

.. math::
   \varepsilon_k \ \equiv \ \|T_k\| \ = \ \beta \|y_k - \tilde{y}_k\|,

corresponding to the local error estimates for three consecutive
steps, :math:`t_{n-3} \to t_{n-2} \to t_{n-1} \to t_n`.  These local
error history values are all initialized to 1 upon program
initialization, to accommodate the few initial time steps of a
calculation where some of these error estimates have not yet been
computed.  With these estimates, ARKODE supports a variety of error
control algorithms, as specified in the subsections below.


.. _ARKODE.Mathematics.Adaptivity.ErrorControl.PID:

PID controller
-----------------

This is the default time adaptivity controller used by the ARKStep and
ERKStep modules.  It derives from those found in :cite:p:`KenCarp:03`, :cite:p:`Sod:98`, :cite:p:`Sod:03` and
:cite:p:`Sod:06`, and uses all three of the local error estimates
:math:`\varepsilon_n`, :math:`\varepsilon_{n-1}` and
:math:`\varepsilon_{n-2}` in determination of a prospective step size,

.. math::
   h' \;=\; h_n\; \varepsilon_n^{-k_1/p}\; \varepsilon_{n-1}^{k_2/p}\;
        \varepsilon_{n-2}^{-k_3/p},

where the constants :math:`k_1`, :math:`k_2` and :math:`k_3` default
to 0.58, 0.21 and 0.1, respectively, and may be modified by the user.
In this estimate, a floor of :math:`\varepsilon > 10^{-10}` is
enforced to avoid division-by-zero errors.



.. _ARKODE.Mathematics.Adaptivity.ErrorControl.PI:

PI controller
----------------------

Like with the previous method, the PI controller derives from those
found in :cite:p:`KenCarp:03`, :cite:p:`Sod:98`, :cite:p:`Sod:03` and :cite:p:`Sod:06`, but it differs in
that it only uses the two most recent step sizes in its adaptivity
algorithm,

.. math::
   h' \;=\; h_n\; \varepsilon_n^{-k_1/p}\; \varepsilon_{n-1}^{k_2/p}.

Here, the default values of :math:`k_1` and :math:`k_2` default
to 0.8 and 0.31, respectively, though they may be changed by the user.



.. _ARKODE.Mathematics.Adaptivity.ErrorControl.I:

I controller
----------------------

This is the standard time adaptivity control algorithm in use by most
publicly-available ODE solver codes.  It bases the prospective time step
estimate entirely off of the current local error estimate,

.. math::
   h' \;=\; h_n\; \varepsilon_n^{-k_1/p}.

By default, :math:`k_1=1`, but that may be modified by the user.




.. _ARKODE.Mathematics.Adaptivity.ErrorControl.eGus:

Explicit Gustafsson controller
---------------------------------

This step adaptivity algorithm was proposed in :cite:p:`Gust:91`, and
is primarily useful with explicit Runge--Kutta methods.
In the notation of our earlier controllers, it has the form

.. math::
   h' \;=\; \begin{cases}
      h_1\; \varepsilon_1^{-1/p}, &\quad\text{on the first step}, \\
      h_n\; \varepsilon_n^{-k_1/p}\;
        \left(\dfrac{\varepsilon_n}{\varepsilon_{n-1}}\right)^{k_2/p}, &
      \quad\text{on subsequent steps}.
   \end{cases}
   :label: ARKODE_expGus

The default values of :math:`k_1` and :math:`k_2` are 0.367 and 0.268,
respectively, and may be modified by the user.




.. _ARKODE.Mathematics.Adaptivity.ErrorControl.iGus:

Implicit Gustafsson controller
---------------------------------

A version of the above controller suitable for implicit Runge--Kutta
methods was introduced in :cite:p:`Gust:94`, and has the form

.. math::
   h' = \begin{cases}
      h_1 \varepsilon_1^{-1/p}, &\quad\text{on the first step}, \\
      h_n \left(\dfrac{h_n}{h_{n-1}}\right) \varepsilon_n^{-k_1/p}
        \left(\dfrac{\varepsilon_n}{\varepsilon_{n-1}}\right)^{-k_2/p}, &
      \quad\text{on subsequent steps}.
   \end{cases}
   :label: ARKODE_impGus

The algorithm parameters default to :math:`k_1 = 0.98` and
:math:`k_2 = 0.95`, but may be modified by the user.




.. _ARKODE.Mathematics.Adaptivity.ErrorControl.ieGus:

ImEx Gustafsson controller
---------------------------------

An ImEx version of these two preceding controllers is also available.
This approach computes the estimates :math:`h'_1` arising from
equation :eq:`ARKODE_expGus` and the estimate :math:`h'_2` arising from
equation :eq:`ARKODE_impGus`, and selects

.. math::
   h' = \frac{h}{|h|}\min\left\{|h'_1|, |h'_2|\right\}.

Here, equation :eq:`ARKODE_expGus` uses :math:`k_1` and
:math:`k_2` with default values of 0.367 and 0.268, while equation
:eq:`ARKODE_impGus` sets both parameters to the input :math:`k_3` that
defaults to 0.95.  All of these values may be modified by the user.



.. _ARKODE.Mathematics.Adaptivity.ErrorControl.User:

User-supplied controller
---------------------------------

Finally, ARKODE's time-stepping modules allow the user to define their
own time step adaptivity function,

.. math::
   h' = H(y, t, h_n, h_{n-1}, h_{n-2}, \varepsilon_n, \varepsilon_{n-1}, \varepsilon_{n-2}, q, p),

to allow for problem-specific choices, or for continued
experimentation with temporal error controllers.





.. _ARKODE.Mathematics.Stability:

Explicit stability
======================

For problems that involve a nonzero explicit component,
i.e. :math:`f^E(t,y) \ne 0` in ARKStep or for any problem in
ERKStep, explicit and ImEx Runge--Kutta methods may benefit from
additional user-supplied information regarding the explicit stability
region.  All ARKODE adaptivity methods utilize estimates of the local
error, and it is often the case that such local error control will be
sufficient for method stability, since unstable steps will typically
exceed the error control tolerances.  However, for problems in which
:math:`f^E(t,y)` includes even moderately stiff components, and
especially for higher-order integration methods, it may occur that
a significant number of attempted steps will exceed the error
tolerances.  While these steps will automatically be recomputed, such
trial-and-error can result in an unreasonable number of failed steps,
increasing the cost of the computation.  In these scenarios, a
stability-based time step controller may also be useful.

Since the maximum stable explicit step for any method depends on the
problem under consideration, in that the value :math:`(h_n\lambda)` must
reside within a bounded stability region, where :math:`\lambda` are
the eigenvalues of the linearized operator :math:`\partial f^E /
\partial y`, information on the maximum stable step size is not
readily available to ARKODE's time-stepping modules.  However, for
many problems such information may be easily obtained through analysis
of the problem itself, e.g. in an advection-diffusion calculation
:math:`f^I` may contain the stiff diffusive components and
:math:`f^E` may contain the comparably nonstiff advection terms.  In
this scenario, an explicitly stable step :math:`h_\text{exp}` would be
predicted as one satisfying the Courant-Friedrichs-Lewy (CFL)
stability condition for the advective portion of the problem,

.. math::
   |h_\text{exp}| < \frac{\Delta x}{|\lambda|}

where :math:`\Delta x` is the spatial mesh size and :math:`\lambda` is
the fastest advective wave speed.

In these scenarios, a user may supply a routine to predict this
maximum explicitly stable step size, :math:`|h_\text{exp}|`.  If a
value for :math:`|h_\text{exp}|` is supplied, it is compared against
the value resulting from the local error controller,
:math:`|h_\text{acc}|`, and the eventual time step used will be
limited accordingly,

.. math::
   h' = \frac{h}{|h|}\min\{c\, |h_\text{exp}|,\, |h_\text{acc}|\}.

Here the explicit stability step factor :math:`c>0` (often called the
"CFL number") defaults to :math:`1/2` but may be modified by the user.




.. _ARKODE.Mathematics.FixedStep:

Fixed time stepping
===================

While both the ARKStep and ERKStep time-stepping modules are designed
for tolerance-based time step adaptivity, they additionally support a
"fixed-step" mode (*note: fixed-step mode is currently required for
the slow time scale in the MRIStep module*).  This mode is typically
used for debugging purposes, for verification against hand-coded
Runge--Kutta methods, or for problems where the time steps should be
chosen based on other problem-specific information.  In this mode,
all internal time step adaptivity is disabled:

* temporal error control is disabled,

* nonlinear or linear solver non-convergence will result in an error
  (instead of a step size adjustment),

* no check against an explicit stability condition is performed.


Additional information on this mode is provided in the sections
:ref:`ARKStep Optional Inputs <ARKODE.Usage.ARKStep.OptionalInputs>`,
:ref:`ERKStep Optional Inputs <ARKODE.Usage.ERKStep.OptionalInputs>`, and
:ref:`MRIStep Optional Inputs <ARKODE.Usage.MRIStep.OptionalInputs>`.





.. _ARKODE.Mathematics.AlgebraicSolvers:

Algebraic solvers
===============================

When solving a problem involving either an implicit component (e.g., in
ARKStep with :math:`f^I(t,y) \ne 0`, or in MRIStep with a solve-decoupled
implicit slow stage), or a non-identity mass matrix (:math:`M(t) \ne I` in
ARKStep), systems of linear or nonlinear algebraic equations must be solved
at each stage and/or step of the method.  This section therefore focuses on
the variety of mathematical methods provided in the ARKODE infrastructure
for such problems, including
:ref:`nonlinear solvers <ARKODE.Mathematics.Nonlinear>`,
:ref:`linear solvers <ARKODE.Mathematics.Linear>`,
:ref:`preconditioners <ARKODE.Mathematics.Preconditioning>`,
:ref:`iterative solver error control <ARKODE.Mathematics.Error>`,
:ref:`implicit predictors <ARKODE.Mathematics.Predictors>`, and techniques
used for simplifying the above solves when using different classes of
:ref:`mass-matrices <ARKODE.Mathematics.MassSolve>`.




.. _ARKODE.Mathematics.Nonlinear:

Nonlinear solver methods
------------------------------------


For the DIRK and ARK methods corresponding to :eq:`ARKODE_IMEX_IVP` and
:eq:`ARKODE_IVP_implicit` in ARKStep, and the solve-decoupled implicit slow
stages :eq:`ARKODE_MRI_delta_c_zero` in MRIStep, an implicit system

.. math::
   G(z_i) = 0
   :label: ARKODE_Residual

must be solved for each implicit stage :math:`z_i`.  In order to
maximize solver efficiency, we define this root-finding problem differently
based on the type of mass-matrix supplied by the user.

* In the case that :math:`M=I` within ARKStep, we define the residual as

  .. math::
     G(z_i) \equiv z_i - h_n A^I_{i,i} f^I(t^I_{n,i}, z_i) - a_i,
     :label: ARKODE_Residual_MeqI

  where we have the data

  .. math::
     a_i \equiv y_{n-1} + h_n \sum_{j=1}^{i-1} \left[
     A^E_{i,j} f^E(t^E_{n,j}, z_j) +
     A^I_{i,j} f^I(t^I_{n,j}, z_j) \right].

* In the case of non-identity mass matrix :math:`M\ne I` within ARKStep, but where
  :math:`M` is independent of :math:`t`, we define the residual as

  .. math::
     G(z_i) \equiv M z_i - h_n A^I_{i,i} f^I(t^I_{n,i}, z_i) - a_i,
     :label: ARKODE_Residual_Mfixed

  where we have the data

  .. math::
     a_i \equiv M y_{n-1} + h_n \sum_{j=1}^{i-1} \left[
     A^E_{i,j} f^E(t^E_{n,j}, z_j) +
     A^I_{i,j} f^I(t^I_{n,j}, z_j) \right].

  .. note::

     This form of residual, as opposed to
     :math:`G(z_i) = z_i - h_n A^I_{i,i} \hat{f}^I(t^I_{n,i}, z_i) - a_i`
     (with :math:`a_i` defined appropriately), removes the need to perform the
     nonlinear solve with right-hand side function :math:`\hat{f}^I=M^{-1}\,f^I`,
     as that would require a linear solve with :math:`M` at *every evaluation* of
     the implicit right-hand side routine.

* In the case of ARKStep with :math:`M` dependent on :math:`t`, we define the residual as

  .. math::
     G(z_i) \equiv M(t^I_{n,i}) (z_i - a_i) - h_n A^I_{i,i} f^I(t^I_{n,i}, z_i)
     :label: ARKODE_Residual_MTimeDep

  where we have the data

  .. math::
     a_i \equiv y_{n-1} + h_n \sum_{j=1}^{i-1} \left[
     A^E_{i,j} \hat{f}^E(t^E_{n,j}, z_j) +
     A^I_{i,j} \hat{f}^I(t^I_{n,j}, z_j) \right].

  .. note::

     As above, this form of the residual is chosen to remove excessive
     mass-matrix solves from the nonlinear solve process.

* Similarly, in MRIStep (that always assumes :math:`M=I`), we have the residual

  .. math::
     G(z_i) \equiv z_i - h^S \left(\sum_{k\geq 0} \frac{\gamma_{i,i}^{\{k\}}}{k+1}\right)
     f^I(t_{n,i}^S, z_i) - a_i = 0
     :label: ARKODE_MRIStep_Residual

  where

  .. math::
     a_i \equiv z_{i-1} + h^S \sum_{j=1}^{i-1} \left(\sum_{k\geq 0}
     \frac{\gamma_{i,j}^{\{k\}}}{k+1}\right)f^I(t_{n,j}^S, z_j).


In each of the above nonlinear residual functions, if :math:`f^I(t,y)` depends
nonlinearly on :math:`y` then :eq:`ARKODE_Residual` corresponds to a nonlinear system
of equations; if instead :math:`f^I(t,y)` depends linearly on :math:`y` then
this is a linear system of equations.

To solve each of the above root-finding problems ARKODE leverages SUNNonlinearSolver
modules from the underlying SUNDIALS infrastructure (see section :numref:`SUNNonlinSol`).
By default, ARKODE selects a variant of :index:`Newton's method`,

.. math::
   z_i^{(m+1)} = z_i^{(m)} + \delta^{(m+1)},
   :label: ARKODE_Newton_iteration

where :math:`m` is the Newton iteration index, and the :index:`Newton
update` :math:`\delta^{(m+1)}` in turn requires the solution of the
:index:`Newton linear system`

.. math::
   {\mathcal A}\left(t^I_{n,i}, z_i^{(m)}\right)\, \delta^{(m+1)} =
   -G\left(z_i^{(m)}\right),
   :label: ARKODE_Newton_system

in which

.. math::
   {\mathcal A}(t,z) \approx M(t) - \gamma J(t,z), \quad
   J(t,z) = \frac{\partial f^I(t,z)}{\partial z}, \quad\text{and}\quad
   \gamma = h_n A^I_{i,i}
   :label: ARKODE_NewtonMatrix

within ARKStep, or

.. math::
   {\mathcal A}(t,z) \approx I - \gamma J(t,z), \quad
   J(t,z) = \frac{\partial f^I(t,z)}{\partial z}, \quad\text{and}\quad
   \gamma = h^S \sum_{k\geq 0} \frac{\gamma_{i,i}^{\{k\}}}{k+1}
   :label: ARKODE_NewtonMatrix_MRIStep

within MRIStep.

In addition to Newton-based nonlinear solvers, the SUNDIALS
SUNNonlinearSolver interface allows solvers of fixed-point type.  These
generally implement a :index:`fixed point iteration` for solving an
implicit stage :math:`z_i`,

.. math::
   z_i^{(m+1)} = \Phi\left(z_i^{(m)}\right) \equiv z_i^{(m)} -
   M(t^I_{n,i})^{-1}\,G\left(z_i^{(m)}\right), \quad m=0,1,\ldots.
   :label: ARKODE_AAFP_iteration

Unlike with Newton-based nonlinear solvers, fixed-point iterations
generally *do not* require the solution of a linear system
involving the Jacobian of :math:`f` at each iteration.

Finally, if the user specifies that :math:`f^I(t,y)` depends linearly on
:math:`y` in ARKStep or MRIStep and if the Newton-based SUNNonlinearSolver
module is used, then the problem :eq:`ARKODE_Residual` will be solved using only a
single Newton iteration.  In this case, an additional user-supplied argument
indicates whether this Jacobian is time-dependent or not, signaling whether the
Jacobian or preconditioner needs to be recomputed at each stage or time step, or
if it can be reused throughout the full simulation.

The optimal choice of solver (Newton vs fixed-point) is highly
problem dependent.  Since fixed-point solvers do not require the
solution of linear systems involving the Jacobian of :math:`f`, each
iteration may be significantly less costly than their Newton
counterparts.  However, this can come at the cost of slower
convergence (or even divergence) in comparison with Newton-like
methods.  While a Newton-based iteration is the default solver in
ARKODE due to its increased robustness on very stiff problems, we
strongly recommend that users also consider the fixed-point solver
when attempting a new problem.

For either the Newton or fixed-point solvers, it is well-known that
both the efficiency and robustness of the algorithm intimately depend
on the choice of a good initial guess.  The initial guess
for these solvers is a prediction :math:`z_i^{(0)}` that is computed
explicitly from previously-computed data (e.g. :math:`y_{n-2}`,
:math:`y_{n-1}`, and :math:`z_j` where :math:`j<i`).  Additional
information on the specific predictor algorithms
is provided in section :numref:`ARKODE.Mathematics.Predictors`.



.. _ARKODE.Mathematics.Linear:

Linear solver methods
------------------------------------

When a Newton-based method is chosen for solving each nonlinear
system, a linear system of equations must be solved at each nonlinear
iteration.  For this solve ARKODE leverages another component of the
shared SUNDIALS infrastructure, the "SUNLinearSolver," described in
section :numref:`SUNLinSol`.   These linear solver modules are grouped
into two categories: matrix-based linear solvers and matrix-free
iterative linear solvers.  ARKODE's interfaces for linear solves of
these types are described in the subsections below.


.. index:: modified Newton iteration

.. _ARKODE.Mathematics.Linear.Direct:

Matrix-based linear solvers
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In the case that a matrix-based linear solver is selected, a *modified
Newton iteration* is utilized.  In a modified Newton iteration, the matrix
:math:`{\mathcal A}` is held fixed for multiple Newton iterations.
More precisely, each Newton iteration is computed from the modified
equation

.. math::
   \tilde{\mathcal A}\left(\tilde{t},\tilde{z}\right)\, \delta^{(m+1)}
   = -G\left(z_i^{(m)}\right),
   :label: ARKODE_modified_Newton_system

in which

.. math::
   \tilde{\mathcal A}(\tilde{t},\tilde{z}) \approx M(\tilde{t}) - \tilde{\gamma} J(\tilde{t},\tilde{z}),
   \quad\text{and}\quad
   \tilde{\gamma} = \tilde{h} A^I_{i,i} \quad\text{(ARKStep)}\\
   :label: ARKODE_modified_NewtonMatrix_ARK

or

.. math::
   \tilde{\mathcal A}(\tilde{t},\tilde{z}) \approx I - \tilde{\gamma} J(\tilde{t},\tilde{z}),
   \quad\text{and}\quad
   \tilde{\gamma} = \tilde{h} \sum_{k\geq 0} \frac{\gamma_{i,i}^{\{k\}}}{k+1}\quad\text{(MRIStep)}.
   :label: ARKODE_modified_NewtonMatrix_MRI

Here, the solution :math:`\tilde{z}`, time :math:`\tilde{t}`, and step
size :math:`\tilde{h}` upon which the modified equation rely, are
merely values of these quantities from a previous iteration.  In other
words, the matrix :math:`\tilde{\mathcal A}` is only computed rarely,
and reused for repeated solves.  As described below in section
:numref:`ARKODE.Mathematics.Linear.Setup`, the frequency at which
:math:`\tilde{\mathcal A}` is recomputed defaults to 20 time steps,
but may be modified by the user.

When using the dense and band SUNMatrix objects for the linear systems
:eq:`ARKODE_modified_Newton_system`, the Jacobian :math:`J` may be supplied
by a user routine, or approximated internally by finite-differences.
In the case of differencing, we use the standard approximation

.. math::
   J_{i,j}(t,z) \approx \frac{f^I_i(t,z+\sigma_j e_j) - f^I_i(t,z)}{\sigma_j},

where :math:`e_j` is the :math:`j`-th unit vector, and the increments
:math:`\sigma_j` are given by

.. math::
   \sigma_j = \max\left\{ \sqrt{U}\, |z_j|, \frac{\sigma_0}{w_j} \right\}.

Here :math:`U` is the unit roundoff, :math:`\sigma_0` is a small
dimensionless value, and :math:`w_j` is the error weight defined in
:eq:`ARKODE_EWT`.  In the dense case, this approach requires :math:`N`
evaluations of :math:`f^I`, one for each column of :math:`J`.  In the
band case, the columns of :math:`J` are computed in groups, using the
Curtis-Powell-Reid algorithm, with the number of :math:`f^I`
evaluations equal to the matrix bandwidth.

We note that with sparse and user-supplied SUNMatrix objects, the
Jacobian *must* be supplied by a user routine.



.. index:: inexact Newton iteration

.. _ARKODE.Mathematics.Linear.Iterative:

Matrix-free iterative linear solvers
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In the case that a matrix-free iterative linear solver is chosen,
an *inexact Newton iteration* is utilized.  Here, the
matrix :math:`{\mathcal A}` is not itself constructed since the
algorithms only require the product of this matrix with a given
vector.  Additionally, each Newton system :eq:`ARKODE_Newton_system` is not
solved completely, since these linear solvers are iterative (hence the
"inexact" in the name). As a result. for these linear solvers
:math:`{\mathcal A}` is applied in a matrix-free manner,

.. math::
   {\mathcal A}(t,z)\, v = M(t)\,v - \gamma\, J(t,z)\, v.

The mass matrix-vector products :math:`Mv` *must* be provided through a
user-supplied routine; the Jacobian matrix-vector products :math:`Jv`
are obtained by either calling an optional user-supplied routine, or
through a finite difference approximation to the directional
derivative:

.. math::
   J(t,z)\,v \approx \frac{f^I(t,z+\sigma v) - f^I(t,z)}{\sigma},

where we use the increment :math:`\sigma = 1/\|v\|` to ensure that
:math:`\|\sigma v\| = 1`.

As with the modified Newton method that reused :math:`{\mathcal A}`
between solves, the inexact Newton iteration may also recompute
the preconditioner :math:`P` infrequently to balance the high costs
of matrix construction and factorization against the reduced
convergence rate that may result from a stale preconditioner.



.. index:: linear solver setup

.. _ARKODE.Mathematics.Linear.Setup:

Updating the linear solver
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In cases where recomputation of the Newton matrix
:math:`\tilde{\mathcal A}` or preconditioner :math:`P` is lagged,
these structures will be recomputed only in the
following circumstances:

* when starting the problem,
* when more than :math:`msbp = 20` steps have been taken since the
  last update (this value may be modified by the user),
* when the value :math:`\tilde{\gamma}` of :math:`\gamma` at the last
  update satisfies :math:`\left|\gamma/\tilde{\gamma} - 1\right| >
  \Delta\gamma_{max} = 0.2` (this value may be modified by the user),
* when a non-fatal convergence failure just occurred,
* when an error test failure just occurred, or
* if the problem is linearly implicit and :math:`\gamma` has
  changed by a factor larger than 100 times machine epsilon.

When an update is forced due to a convergence failure, an update of
:math:`\tilde{\mathcal A}` or :math:`P` may or may not involve a
re-evaluation of :math:`J` (in :math:`\tilde{\mathcal A}`) or of
Jacobian data (in :math:`P`), depending on whether errors in the
Jacobian were the likely cause of the failure.  More generally, the
decision is made to re-evaluate :math:`J` (or instruct the user to
update :math:`P`) when:

* starting the problem,
* more than :math:`msbj=50` steps have been taken since the last evaluation,
* a convergence failure occurred with an outdated matrix, and the
  value :math:`\tilde{\gamma}` of :math:`\gamma` at the last update
  satisfies :math:`\left|\gamma/\tilde{\gamma} - 1\right| > 0.2`,
* a convergence failure occurred that forced a step size reduction, or
* if the problem is linearly implicit and :math:`\gamma` has
  changed by a factor larger than 100 times machine epsilon.


However, for linear solvers and preconditioners that do not
rely on costly matrix construction and factorization operations
(e.g. when using a geometric multigrid method as preconditioner), it
may be more efficient to update these structures more frequently than
the above heuristics specify, since the increased rate of
linear/nonlinear solver convergence may more than account for the
additional cost of Jacobian/preconditioner construction.  To this end,
a user may specify that the system matrix :math:`{\mathcal A}` and/or
preconditioner :math:`P` should be recomputed more frequently.

As will be further discussed in section :numref:`ARKODE.Mathematics.Preconditioning`,
in the case of most Krylov methods, preconditioning may be applied on the
left, right, or on both sides of :math:`{\mathcal A}`, with user-supplied
routines for the preconditioner setup and solve operations.




.. _ARKODE.Mathematics.Error:

Iteration Error Control
------------------------------------


.. _ARKODE.Mathematics.Error.Nonlinear:

Nonlinear iteration error control
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

ARKODE provides a customized stopping test to the SUNNonlinearSolver
module used for solving equation :eq:`ARKODE_Residual`.  This test is related
to the temporal local error test, with the goal of keeping the
nonlinear iteration errors from interfering with local error control.
Denoting the final computed value of each stage solution as
:math:`z_i^{(m)}`, and the true stage solution solving :eq:`ARKODE_Residual`
as :math:`z_i`, we want to ensure that the iteration error
:math:`z_i - z_i^{(m)}` is "small" (recall that a norm less than 1 is
already considered within an acceptable tolerance).

To this end, we first estimate the linear convergence rate :math:`R_i`
of the nonlinear iteration.  We initialize :math:`R_i=1`, and reset it
to this value whenever :math:`\tilde{\mathcal A}` or :math:`P` are
updated.  After computing a nonlinear correction :math:`\delta^{(m)} =
z_i^{(m)} - z_i^{(m-1)}`, if :math:`m>0` we update :math:`R_i` as

.. math::
   R_i \leftarrow \max\left\{ c_r R_i, \left\|\delta^{(m)}\right\| / \left\|\delta^{(m-1)}\right\| \right\}.
   :label: ARKODE_NonlinearCRate

where the default factor :math:`c_r=0.3` is user-modifiable.

Let :math:`y_n^{(m)}` denote the time-evolved solution constructed
using our approximate nonlinear stage solutions, :math:`z_i^{(m)}`,
and let :math:`y_n^{(\infty)}` denote the time-evolved solution
constructed using *exact* nonlinear stage solutions.  We then use the
estimate

.. math::
   \left\| y_n^{(\infty)} - y_n^{(m)} \right\| \approx
   \max_i \left\| z_i^{(m+1)} - z_i^{(m)} \right\| \approx
   \max_i R_i \left\| z_i^{(m)} - z_i^{(m-1)} \right\| =
   \max_i R_i \left\| \delta^{(m)} \right\|.

Therefore our convergence (stopping) test for the nonlinear iteration
for each stage is

.. math::
   R_i \left\|\delta^{(m)} \right\| < \epsilon,
   :label: ARKODE_NonlinearTolerance

where the factor :math:`\epsilon` has default value 0.1.  We default
to a maximum of 3 nonlinear iterations.  We also declare the
nonlinear iteration to be divergent if any of the ratios

.. math::
   `\|\delta^{(m)}\| / \|\delta^{(m-1)}\| > r_{div}`
   :label: ARKODE_NonlinearDivergence

with :math:`m>0`, where :math:`r_{div}` defaults to 2.3.
If convergence fails in the nonlinear solver with :math:`{\mathcal A}`
current (i.e., not lagged), we reduce the step size :math:`h_n` by a
factor of :math:`\eta_{cf}=0.25`.  The integration will be halted after
:math:`max_{ncf}=10` convergence failures, or if a convergence failure
occurs with :math:`h_n = h_\text{min}`.  However, since the nonlinearity
of :eq:`ARKODE_Residual` may vary significantly based on the problem under
consideration, these default constants may all be modified by the user.



.. _ARKODE.Mathematics.Error.Linear:

Linear iteration error control
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

When a Krylov method is used to solve the linear Newton systems
:eq:`ARKODE_Newton_system`, its errors must also be controlled.  To this end,
we approximate the linear iteration error in the solution vector
:math:`\delta^{(m)}` using the preconditioned residual vector,
e.g. :math:`r = P{\mathcal A}\delta^{(m)} + PG` for the case of left
preconditioning (the role of the preconditioner is further elaborated
in the next section).  In an attempt to ensure that the linear
iteration errors do not interfere with the nonlinear solution error
and local time integration error controls, we require that the norm of
the preconditioned linear residual satisfies

.. math::
   \|r\| \le \frac{\epsilon_L \epsilon}{10}.
   :label: ARKODE_LinearTolerance

Here :math:`\epsilon` is the same value as that is used above for the
nonlinear error control.  The factor of 10 is used to ensure that the
linear solver error does not adversely affect the nonlinear solver
convergence.  Smaller values for the parameter :math:`\epsilon_L` are
typically useful for strongly nonlinear or very stiff ODE systems,
while easier ODE systems may benefit from a value closer to 1.  The
default value is :math:`\epsilon_L = 0.05`, which may be modified by
the user.  We note that for linearly
implicit problems the tolerance :eq:`ARKODE_LinearTolerance` is similarly
used for the single Newton iteration.




.. _ARKODE.Mathematics.Preconditioning:

Preconditioning
------------------------------------

When using an inexact Newton method to solve the nonlinear system
:eq:`ARKODE_Residual`, an iterative method is used repeatedly to solve
linear systems of the form :math:`{\mathcal A}x = b`, where :math:`x` is a
correction vector and :math:`b` is a residual vector.  If this
iterative method is one of the scaled preconditioned iterative linear
solvers supplied with SUNDIALS, their efficiency may benefit
tremendously from preconditioning. A system :math:`{\mathcal A}x=b`
can be preconditioned using any one of:

.. math::
   (P^{-1}{\mathcal A})x = P^{-1}b & \qquad\text{[left preconditioning]}, \\
   ({\mathcal A}P^{-1})Px = b  & \qquad\text{[right preconditioning]}, \\
   (P_L^{-1} {\mathcal A} P_R^{-1}) P_R x = P_L^{-1}b & \qquad\text{[left and right
   preconditioning]}.

These Krylov iterative methods are then applied to a system with the
matrix :math:`P^{-1}{\mathcal A}`, :math:`{\mathcal A}P^{-1}`, or
:math:`P_L^{-1} {\mathcal A} P_R^{-1}`, instead of :math:`{\mathcal
A}`.  In order to improve the convergence of the Krylov iteration, the
preconditioner matrix :math:`P`, or the product :math:`P_L P_R` in the
third case, should in some sense approximate the system matrix
:math:`{\mathcal A}`.  Simultaneously, in order to be
cost-effective the matrix :math:`P` (or matrices :math:`P_L` and
:math:`P_R`) should be reasonably efficient to evaluate and solve.
Finding an optimal point in this trade-off between rapid
convergence and low cost can be quite challenging.  Good choices are
often problem-dependent (for example, see :cite:p:`BrHi:89` for an
extensive study of preconditioners for reaction-transport systems).

Most of the iterative linear solvers supplied with SUNDIALS allow for
all three types of preconditioning (left, right or both), although for
non-symmetric matrices :math:`{\mathcal A}` we know of few situations
where preconditioning on both sides is superior to preconditioning on
one side only (with the product :math:`P = P_L P_R`).  Moreover, for a
given preconditioner matrix, the merits of left vs. right
preconditioning are unclear in general, so we recommend that the user
experiment with both choices.  Performance can differ between these
since the inverse of the left preconditioner is included in the linear
system residual whose norm is being tested in the Krylov algorithm.
As a rule, however, if the preconditioner is the product of two
matrices, we recommend that preconditioning be done either on the left
only or the right only, rather than using one factor on each
side.  An exception to this rule is the PCG solver, that itself
assumes a symmetric matrix :math:`{\mathcal A}`, since the PCG
algorithm in fact applies the single preconditioner matrix :math:`P`
in both left/right fashion as :math:`P^{-1/2} {\mathcal A} P^{-1/2}`.

Typical preconditioners are based on approximations
to the system Jacobian, :math:`J = \partial f^I / \partial y`.  Since
the Newton iteration matrix involved is :math:`{\mathcal A} = M -
\gamma J`, any approximation :math:`\bar{J}` to :math:`J` yields a
matrix that is of potential use as a preconditioner, namely :math:`P =
M - \gamma \bar{J}`. Because the Krylov iteration occurs within a
Newton iteration and further also within a time integration, and since
each of these iterations has its own test for convergence, the
preconditioner may use a very crude approximation, as long as it
captures the dominant numerical features of the system.  We have
found that the combination of a preconditioner with the Newton-Krylov
iteration, using even a relatively poor approximation to the Jacobian,
can be surprisingly superior to using the same matrix without Krylov
acceleration (i.e., a modified Newton iteration), as well as to using
the Newton-Krylov method with no preconditioning.




.. _ARKODE.Mathematics.Predictors:

Implicit predictors
------------------------------------

For problems with implicit components, a prediction algorithm is
employed for constructing the initial guesses for each implicit
Runge--Kutta stage, :math:`z_i^{(0)}`.  As is well-known with nonlinear
solvers, the selection of a good initial guess can have dramatic
effects on both the speed and robustness of the solve, making the
difference between rapid quadratic convergence versus divergence of
the iteration.  To this end, a variety of prediction algorithms are
provided.  In each case, the stage guesses :math:`z_i^{(0)}` are
constructed explicitly using readily-available information, including
the previous step solutions :math:`y_{n-1}` and :math:`y_{n-2}`, as
well as any previous stage solutions :math:`z_j, \quad j<i`.  In most
cases, prediction is performed by constructing an interpolating
polynomial through existing data, which is then evaluated at the
desired stage time to provide an inexpensive but (hopefully)
reasonable prediction of the stage solution.  Specifically, for most
Runge--Kutta methods each stage solution satisfies

.. math::
   z_i \approx y(t^I_{n,i}),

(similarly for MRI methods :math:`z_i \approx y(t^S_{n,i})`),
so by constructing an interpolating polynomial :math:`p_q(t)` through
a set of existing data, the initial guess at stage solutions may be
approximated as

.. math::
   z_i^{(0)} = p_q(t^I_{n,i}).
   :label: ARKODE_extrapolant

As the stage times for MRI stages and implicit ARK and DIRK stages usually
have non-negative abscissae (i.e., :math:`c_j^I > 0`), it is typically the
case that :math:`t^I_{n,j}` (resp., :math:`t^S_{n,j}`) is outside of the
time interval containing the data used to construct :math:`p_q(t)`, hence
:eq:`ARKODE_extrapolant` will correspond to an extrapolant instead of an
interpolant.  The dangers of using a polynomial interpolant to extrapolate
values outside the interpolation interval are well-known, with higher-order
polynomials and predictions further outside the interval resulting in the
greatest potential inaccuracies.

The prediction algorithms available in ARKODE therefore
construct a variety of interpolants :math:`p_q(t)`, having
different polynomial order and using different interpolation data, to
support "optimal" choices for different types of problems, as
described below.  We note that due to the structural similarities between
implicit ARK and DIRK stages in ARKStep, and solve-decoupled implicit stages
in MRIStep, we use the ARKStep notation throughout the remainder of this
section, but each statement equally applies to MRIStep (unless otherwise noted).


.. _ARKODE.Mathematics.Predictors.Trivial:

Trivial predictor
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The so-called "trivial predictor" is given by the formula

.. math::

   p_0(t) = y_{n-1}.

While this piecewise-constant interpolant is clearly not a highly
accurate candidate for problems with time-varying solutions, it is
often the most robust approach for highly stiff problems, or for
problems with implicit constraints whose violation may cause illegal
solution values (e.g. a negative density or temperature).


.. _ARKODE.Mathematics.Predictors.Max:

Maximum order predictor
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

At the opposite end of the spectrum, ARKODE's interpolation modules
discussed in section :numref:`ARKODE.Mathematics.Interpolation`
can be used to construct a higher-order polynomial interpolant, :math:`p_q(t)`.
The implicit stage predictor is computed through evaluating the
highest-degree-available interpolant at each stage time :math:`t^I_{n,i}`.



.. _ARKODE.Mathematics.Predictors.Decreasing:

Variable order predictor
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This predictor attempts to use higher-degree polynomials
:math:`p_q(t)` for predicting earlier stages, and lower-degree
interpolants for later stages.  It uses the same interpolation module
as described above, but chooses the polynomial degree adaptively based on the
stage index :math:`i`, under the assumption that the
stage times are increasing, i.e. :math:`c^I_j < c^I_k` for
:math:`j<k`:

.. math::
   q_i = \max\{ q_\text{max} - i + 1,\; 1 \}, \quad i=1,\ldots,s.



.. _ARKODE.Mathematics.Predictors.Cutoff:

Cutoff order predictor
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This predictor follows a similar idea as the previous algorithm, but
monitors the actual stage times to determine the polynomial
interpolant to use for prediction.  Denoting :math:`\tau = c_i^I
\dfrac{h_n}{h_{n-1}}`, the polynomial degree :math:`q_i` is chosen as:

.. math::
   q_i = \begin{cases}
      q_\text{max}, & \text{if}\quad \tau < \tfrac12,\\
      1, & \text{otherwise}.
   \end{cases}



.. _ARKODE.Mathematics.Predictors.Bootstrap:

Bootstrap predictor (:math:`M=I` only) -- **deprecated**
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This predictor does not use any information from the preceding
step, instead using information only within the current step
:math:`[t_{n-1},t_n]`.  In addition to using the solution and ODE
right-hand side function, :math:`y_{n-1}` and
:math:`f(t_{n-1},y_{n-1})`, this approach uses the right-hand
side from a previously computed stage solution in the same step,
:math:`f(t_{n-1}+c^I_j h,z_j)` to construct a quadratic Hermite
interpolant for the prediction.  If we define the constants
:math:`\tilde{h} = c^I_j h` and :math:`\tau = c^I_i h`, the predictor
is given by

.. math::

   z_i^{(0)} = y_{n-1} + \left(\tau - \frac{\tau^2}{2\tilde{h}}\right)
      f(t_{n-1},y_{n-1}) + \frac{\tau^2}{2\tilde{h}} f(t_{n-1}+\tilde{h},z_j).

For stages without a nonzero preceding stage time,
i.e. :math:`c^I_j\ne 0` for :math:`j<i`, this method reduces to using
the trivial predictor :math:`z_i^{(0)} = y_{n-1}`.  For stages having
multiple preceding nonzero :math:`c^I_j`, we choose the stage having
largest :math:`c^I_j` value, to minimize the level of extrapolation
used in the prediction.

We note that in general, each stage solution :math:`z_j` has
significantly worse accuracy than the time step solutions
:math:`y_{n-1}`, due to the difference between the *stage order* and
the *method order* in Runge--Kutta methods.  As a result, the accuracy
of this predictor will generally be rather limited, but it is
provided for problems in which this increased stage error is better
than the effects of extrapolation far outside of the previous time
step interval :math:`[t_{n-2},t_{n-1}]`.

Although this approach could be used with non-identity mass matrix, support for
that mode is not currently implemented, so selection of this predictor in the
case of a non-identity mass matrix will result in use of the trivial predictor.

.. note::
   This predictor has been deprecated, and will be removed from a future release.



.. _ARKODE.Mathematics.Predictors.MinimumCorrection:

Minimum correction predictor (ARKStep, :math:`M=I` only) -- **deprecated**
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The final predictor is not interpolation based; instead it
utilizes all existing stage information from the current step to
create a predictor containing all but the current stage solution.
Specifically, as discussed in equations :eq:`ARKODE_ARK` and :eq:`ARKODE_Residual`,
each stage solves a nonlinear equation

.. math::
   z_i &= y_{n-1} + h_n \sum_{j=1}^{i-1} A^E_{i,j} f^E(t^E_{n,j}, z_j)
   + h_n \sum_{j=1}^{i}   A^I_{i,j} f^I(t^I_{n,j}, z_j), \\
   \Leftrightarrow \qquad \qquad & \\
   G(z_i) &\equiv z_i - h_n A^I_{i,i} f^I(t^I_{n,i}, z_i) - a_i = 0.

This prediction method merely computes the predictor :math:`z_i` as

.. math::
   z_i &= y_{n-1} + h_n \sum_{j=1}^{i-1} A^E_{i,j} f^E(t^E_{n,j}, z_j)
                 + h_n \sum_{j=1}^{i-1}  A^I_{i,j} f^I(t^I_{n,j}, z_j), \\
   \Leftrightarrow \quad \qquad & \\
   z_i &= a_i.

Again, although this approach could be used with non-identity mass matrix, support
for that mode is not currently implemented, so selection of this predictor in the
case of a non-identity mass matrix will result in use of the trivial predictor.

.. note::
   This predictor has been deprecated, and will be removed from a future release.





.. _ARKODE.Mathematics.MassSolve:

Mass matrix solver (ARKStep only)
------------------------------------

Within the ARKStep algorithms described above, there are multiple
locations where a matrix-vector product

.. math::
   b = M v
   :label: ARKODE_mass_multiply

or a linear solve

.. math::
   x = M^{-1} b
   :label: ARKODE_mass_solve

is required.

Of course, for problems in which :math:`M=I` both of these operators
are trivial.  However for problems with non-identity mass matrix,
these linear solves :eq:`ARKODE_mass_solve` may be handled using
any valid SUNLinearSolver module, in the same manner as described in the
section :numref:`ARKODE.Mathematics.Linear` for solving the linear Newton
systems.

For ERK methods involving non-identity mass matrix, even though
calculation of individual stages does not require an algebraic solve,
both of the above operations (matrix-vector product, and mass matrix
solve) may be required within each time step.  Therefore, for these
users we recommend reading the rest of this section as it pertains to
ARK methods, with the obvious simplification that since :math:`f^E=f`
and :math:`f^I=0` no Newton or fixed-point nonlinear solve, and no
overall system linear solve, is involved in the solution process.

At present, for DIRK and ARK problems using a matrix-based solver for
the Newton nonlinear iterations, the type of matrix (dense, band,
sparse, or custom) for the Jacobian matrix :math:`J` must match the
type of mass matrix :math:`M`, since these are combined to form the
Newton system matrix :math:`\tilde{\mathcal A}`.  When matrix-based
methods are employed, the user must supply a routine to compute
:math:`M(t)` in the appropriate form to match the structure of
:math:`{\mathcal A}`, with a user-supplied routine of type
:c:func:`ARKLsMassFn()`.  This matrix structure is used internally to
perform any requisite mass matrix-vector products :eq:`ARKODE_mass_multiply`.

When matrix-free methods are selected, a routine must be supplied to
perform the mass-matrix-vector product, :math:`Mv`.  As with iterative
solvers for the Newton systems, preconditioning may be applied to aid
in solution of the mass matrix systems :eq:`ARKODE_mass_solve`.  When using an
iterative mass matrix linear solver, we require that the norm of the
preconditioned linear residual satisfies

.. math::
   \|r\| \le \epsilon_L \epsilon,
   :label: ARKODE_MassLinearTolerance

where again, :math:`\epsilon` is the nonlinear solver tolerance
parameter from :eq:`ARKODE_NonlinearTolerance`.  When using iterative system
and mass matrix linear solvers, :math:`\epsilon_L` may be specified
separately for both tolerances :eq:`ARKODE_LinearTolerance` and
:eq:`ARKODE_MassLinearTolerance`.


In the algorithmic descriptions above there are five locations
where a linear solve of the form :eq:`ARKODE_mass_solve` is required: (a) at each
iteration of a fixed-point nonlinear solve, (b) in computing the
Runge--Kutta right-hand side vectors :math:`\hat{f}^E_i` and
:math:`\hat{f}^I_i`, (c) in constructing the time-evolved solution
:math:`y_n`, (d) in estimating the local temporal truncation error, and (e)
in constructing predictors for the implicit solver iteration (see section
:numref:`ARKODE.Mathematics.Predictors.Max`).  We note that different nonlinear
solver approaches (i.e., Newton vs fixed-point) and different types of
mass matrices (i.e., time-dependent versus fixed) result in different
subsets of the above operations.  We discuss each of these in the bullets below.

* When using a fixed-point nonlinear solver, at each fixed-point iteration
  we must solve

  .. math::
     M(t^I_{n,i})\, z_i^{(m+1)} = G\left(z_i^{(m)}\right), \quad m=0,1,\ldots

  for the new fixed-point iterate, :math:`z_i^{(m+1)}`.

* In the case of a time-dependent mass matrix, to construct the Runge--Kutta
  right-hand side vectors we must solve

  .. math::
     M(t^E_{n,i}) \hat{f}^{E}_i \ = \ f^{E}(t^E_{n,i},z_i)
     \quad\text{and}\quad
     M(t^I_{n,i}) \hat{f}^{I}_j \ = \ f^{I}(t^I_{n,i},z_i)

  for the vectors :math:`\hat{f}^{E}_i` and :math:`\hat{f}^{I}_i`.

* For fixed mass matrices, we construct the time-evolved solution :math:`y_n`
  from equation :eq:`ARKODE_ARK` by solving

  .. math::
     &M y_n \ = \ M y_{n-1} + h_n \sum_{i=1}^{s} \left( b^E_i f^E(t^E_{n,i}, z_i)
                   + b^I_i f^I(t^I_{n,i}, z_i)\right), \\
     \Leftrightarrow \qquad & \\
     &M (y_n -y_{n-1}) \ = \ h_n \sum_{i=1}^{s} \left(b^E_i f^E(t^E_{n,i}, z_i)
                   + b^I_i f^I(t^I_{n,i}, z_i)\right), \\
     \Leftrightarrow \qquad & \\
     &M \nu \ = \ h_n \sum_{i=1}^{s} \left(b^E_i f^E(t^E_{n,i}, z_i)
                   + b^I_i f^I(t^I_{n,i}, z_i)\right),

  for the update :math:`\nu = y_n - y_{n-1}`.

  Similarly, we compute the local temporal error
  estimate :math:`T_n` from equation :eq:`ARKODE_LTE` by solving systems of the form

  .. math::
     M\, T_n = h \sum_{i=1}^{s} \left[
     \left(b^E_i - \tilde{b}^E_i\right) f^E(t^E_{n,i}, z_i) +
     \left(b^I_i - \tilde{b}^I_i\right) f^I(t^I_{n,i}, z_i) \right].
     :label: ARKODE_mass_solve_LTE

* For problems with either form of non-identity mass matrix, in constructing
  dense output and implicit predictors of degree 2 or higher (see the
  section :numref:`ARKODE.Mathematics.Predictors.Max` above), we compute the derivative
  information :math:`\hat{f}_k` from the equation

  .. math::
     M(t_n) \hat{f}_n = f^E(t_n, y_n) + f^I(t_n, y_n).

In total, for problems with fixed mass matrix, we require only
two mass-matrix linear solves :eq:`ARKODE_mass_solve` per attempted time step,
with one more upon completion of a time step that meets the solution accuracy
requirements.  When fixed time-stepping is used (:math:`h_n=h`), the
solve :eq:`ARKODE_mass_solve_LTE` is not performed at each attempted step.

Similarly, for problems with time-dependent mass matrix, we require
:math:`2s` mass-matrix linear solves :eq:`ARKODE_mass_solve` per attempted step,
where :math:`s` is the number of stages in the ARK method (only half of
these are required for purely explicit or purely implicit problems,
:eq:`ARKODE_IVP_explicit` or :eq:`ARKODE_IVP_implicit`), with one more upon completion of
a time step that meets the solution accuracy requirements.

In addition to the above totals, when using a fixed-point nonlinear solver
(assumed to require :math:`m` iterations), we will need an additional
:math:`ms` mass-matrix linear solves :eq:`ARKODE_mass_solve` per attempted time
step (but zero linear solves with the system Jacobian).



.. _ARKODE.Mathematics.Rootfinding:

Rootfinding
===============

All of the time-stepping modules in ARKODE also support a rootfinding
feature.  This means that, while integrating the IVP :eq:`ARKODE_IVP`, these
can also find the roots of a set of user-defined functions
:math:`g_i(t,y)` that depend on :math:`t` and the solution vector
:math:`y = y(t)`. The number of these root functions is arbitrary, and
if more than one :math:`g_i` is found to have a root in any given
interval, the various root locations are found and reported in the
order that they occur on the :math:`t` axis, in the direction of
integration.

Generally, this rootfinding feature finds only roots of odd
multiplicity, corresponding to changes in sign of :math:`g_i(t,
y(t))`, denoted :math:`g_i(t)` for short. If a user root function has
a root of even multiplicity (no sign change), it will almost certainly
be missed due to the realities of floating-point arithmetic.  If such
a root is desired, the user should reformulate the root function so
that it changes sign at the desired root.

The basic scheme used is to check for sign changes of any
:math:`g_i(t)` over each time step taken, and then (when a sign change
is found) to home in on the root (or roots) with a modified secant
method :cite:p:`HeSh:80`.  In addition, each time :math:`g` is
evaluated, ARKODE checks to see if :math:`g_i(t) = 0` exactly, and if
so it reports this as a root.  However, if an exact zero of any
:math:`g_i` is found at a point :math:`t`, ARKODE computes
:math:`g(t+\delta)` for a small increment :math:`\delta`, slightly
further in the direction of integration, and if any
:math:`g_i(t+\delta) = 0` also, ARKODE stops and reports an
error. This way, each time ARKODE takes a time step, it is guaranteed
that the values of all :math:`g_i` are nonzero at some past value of
:math:`t`, beyond which a search for roots is to be done.

At any given time in the course of the time-stepping, after suitable
checking and adjusting has been done, ARKODE has an interval
:math:`(t_\text{lo}, t_\text{hi}]` in which roots of the
:math:`g_i(t)` are to be sought, such that :math:`t_\text{hi}` is
further ahead in the direction of integration, and all
:math:`g_i(t_\text{lo}) \ne 0`.  The endpoint :math:`t_\text{hi}` is
either :math:`t_n`, the end of the time step last taken, or the next
requested output time :math:`t_\text{out}` if this comes sooner. The
endpoint :math:`t_\text{lo}` is either :math:`t_{n-1}`, or the last
output time :math:`t_\text{out}` (if this occurred within the last
step), or the last root location (if a root was just located within
this step), possibly adjusted slightly toward :math:`t_n` if an exact
zero was found. The algorithm checks :math:`g(t_\text{hi})` for zeros, and
it checks for sign changes in :math:`(t_\text{lo}, t_\text{hi})`. If no sign
changes are found, then either a root is reported (if some
:math:`g_i(t_\text{hi}) = 0`) or we proceed to the next time interval
(starting at :math:`t_\text{hi}`). If one or more sign changes were found,
then a loop is entered to locate the root to within a rather tight
tolerance, given by

.. math::
   \tau = 100\, U\, (|t_n| + |h|)\qquad (\text{where}\; U = \text{unit roundoff}).

Whenever sign changes are seen in two or more root functions, the one
deemed most likely to have its root occur first is the one with the
largest value of
:math:`\left|g_i(t_\text{hi})\right| / \left| g_i(t_\text{hi}) - g_i(t_\text{lo})\right|`,
corresponding to the closest to :math:`t_\text{lo}` of the secant method
values. At each pass through the loop, a new value :math:`t_\text{mid}` is
set, strictly within the search interval, and the values of
:math:`g_i(t_\text{mid})` are checked. Then either :math:`t_\text{lo}` or
:math:`t_\text{hi}` is reset to :math:`t_\text{mid}` according to which
subinterval is found to have the sign change. If there is none in
:math:`(t_\text{lo}, t_\text{mid})` but some :math:`g_i(t_\text{mid}) = 0`, then that
root is reported. The loop continues until :math:`\left|t_\text{hi} -
t_\text{lo} \right| < \tau`, and then the reported root location is
:math:`t_\text{hi}`.  In the loop to locate the root of :math:`g_i(t)`, the
formula for :math:`t_\text{mid}` is

.. math::
   t_\text{mid} = t_\text{hi} -
   \frac{g_i(t_\text{hi}) (t_\text{hi} - t_\text{lo})}{g_i(t_\text{hi}) - \alpha g_i(t_\text{lo})} ,

where :math:`\alpha` is a weight parameter. On the first two passes
through the loop, :math:`\alpha` is set to 1, making :math:`t_\text{mid}`
the secant method value. Thereafter, :math:`\alpha` is reset according
to the side of the subinterval (low vs high, i.e. toward
:math:`t_\text{lo}` vs toward :math:`t_\text{hi}`) in which the sign change was
found in the previous two passes. If the two sides were opposite,
:math:`\alpha` is set to 1. If the two sides were the same, :math:`\alpha`
is halved (if on the low side) or doubled (if on the high side). The
value of :math:`t_\text{mid}` is closer to :math:`t_\text{lo}` when
:math:`\alpha < 1` and closer to :math:`t_\text{hi}` when :math:`\alpha > 1`.
If the above value of :math:`t_\text{mid}` is within :math:`\tau /2` of
:math:`t_\text{lo}` or :math:`t_\text{hi}`, it is adjusted inward, such that its
fractional distance from the endpoint (relative to the interval size)
is between 0.1 and 0.5 (with 0.5 being the midpoint), and the actual
distance from the endpoint is at least :math:`\tau/2`.

Finally, we note that when running in parallel, ARKODE's rootfinding
module assumes that the entire set of root defining functions
:math:`g_i(t,y)` is replicated on every MPI rank.  Since in these
cases the vector :math:`y` is distributed across ranks, it is the
user's responsibility to perform any necessary communication to ensure
that :math:`g_i(t,y)` is identical on each rank.


.. _ARKODE.Mathematics.InequalityConstraints:

Inequality Constraints
=======================

The ARKStep and ERKStep modules in ARKODE permit the user to impose optional
inequality constraints on individual components of the solution vector :math:`y`.
Any of the following four constraints can be imposed: :math:`y_i > 0`, :math:`y_i < 0`,
:math:`y_i \geq 0`, or :math:`y_i \leq 0`. The constraint satisfaction is tested
after a successful step and before the error test. If any constraint fails, the
step size is reduced and a flag is set to update the Jacobian or preconditioner
if applicable. Rather than cutting the step size by some arbitrary factor,
ARKODE estimates a new step size :math:`h'` using a linear approximation of the
components in :math:`y` that failed the constraint test (including a safety
factor of 0.9 to cover the strict inequality case). If a step fails to satisfy
the constraints 10 times (a value which may be modified by the user) within a
step attempt, or fails with the minimum step size, then the integration is halted
and an error is returned. In this case the user may need to employ other
strategies as discussed in :numref:`ARKODE.Usage.ARKStep.Tolerances` and
:numref:`ARKODE.Usage.ERKStep.Tolerances` to satisfy the inequality constraints.
