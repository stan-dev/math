.. ----------------------------------------------------------------
   SUNDIALS Copyright Start
   Copyright (c) 2002-2021, Lawrence Livermore National Security
   and Southern Methodist University.
   All rights reserved.

   See the top-level LICENSE and NOTICE files for details.

   SPDX-License-Identifier: BSD-3-Clause
   SUNDIALS Copyright End
   ----------------------------------------------------------------

.. _KINSOL.Mathematics:

***************************
Mathematical Considerations
***************************

KINSOL solves nonlinear algebraic systems in real :math:`N`-space.

Using Newton’s method, or the Picard iteration, one can solve

.. math::
  F(u) = 0 \, , \quad F:\mathbb{R}^N \rightarrow \mathbb{R}^N \, ,
  :label: KIN_nonlinear_system

given an initial guess :math:`u_0`. Using a fixed-point iteration, the convergence of which can be improved with
Anderson acceleration, one can solve

.. math::
  G(u) = u \, , \quad G:\mathbb{R}^N \rightarrow \mathbb{R}^N \, ,
  :label: KIN_fixed-point_system

given an initial guess :math:`u_0`.

Basic Newton iteration
----------------------

Depending on the linear solver used, KINSOL can employ either an Inexact Newton
method :cite:p:`Bro:87,BrSa:90,DES:82,DeSc:96,Kel:95`, or a Modified Newton
method. At the highest level, KINSOL implements the following iteration scheme:

#. Set :math:`u_0 =` an initial guess

#. For :math:`n = 0, 1, 2,...` until convergence do:

   .. _Newton:

   a. Solve :math:`J(u_n)\delta_n = -F(u_n)`

   b. Set :math:`u_{n+1} = u_n + \lambda \delta_n`, :math:`0 < \lambda \leq 1`

   c. Test for convergence

Here, :math:`u_n` is the :math:`n`\ th iterate to :math:`u`, and :math:`J(u) = F'(u)` is the system Jacobian. At each
stage in the iteration process, a scalar multiple of the step :math:`\delta_n`, is added to :math:`u_n` to produce a new
iterate, :math:`u_{n+1}`. A test for convergence is made before the iteration continues.

Newton method variants
----------------------

For solving the linear system given in step (:ref:`2a <Newton>`), KINSOL provides several choices,
including the option of a user-supplied linear solver module. The linear solver modules distributed with SUNDIALS
are organized in two families, a *direct* family comprising direct linear solvers for dense, banded, or sparse matrices
and a *spils* family comprising scaled preconditioned iterative (Krylov) linear solvers. The methods offered through
these modules are as follows:

-  dense direct solvers, using either an internal implementation or a BLAS/LAPACK implementation (serial
   or threaded vector modules only),

-  band direct solvers, using either an internal implementation or a BLAS/LAPACK implementation (serial
   or threaded vector modules only),

-  sparse direct solver interfaces to various libraries, including KLU :cite:p:`DaPa:10, KLU_site`,
   SuperLU_MT :cite:p:`Li:05,DGL:99,SuperLUMT_site`, SuperLU_Dist
   :cite:p:`GDL:07,LD:03,SLUUG:99,SuperLUDIST_site`, and cuSPARSE :cite:p:`cuSPARSE_site` [Note that users
   will need to download and install the relevant external packages independent of KINSOL],

-  SPGMR, a scaled preconditioned GMRES (Generalized Minimal Residual method) solver,

-  SPFGMR, a scaled preconditioned FGMRES (Flexible Generalized Minimal Residual method) solver,

-  SPBCG, a scaled preconditioned Bi-CGStab (Bi-Conjugate Gradient Stable method) solver,

-  SPTFQMR, a scaled preconditioned TFQMR (Transpose-Free Quasi-Minimal Residual method) solver, or

-  PCG, a scaled preconditioned CG (Conjugate Gradient method) solver.

When using a direct linear solver, the linear system in :ref:`2a <Newton>` is
solved exactly, thus resulting in a Modified Newton method (the Jacobian matrix
is normally out of date; see below). Note that KINSOL allows the user to enforce
a Jacobian evaluation at each iteration thus allowing for an Exact Newton
iteration. Note that each direct linear solver is only compatible with a subset of
vector representations (see :numref:`SUNLinSol.API.Compatibility` for details).

When using an iterative linear solver, the linear system in (:ref:`2a <Newton>`)
is solved only approximately, thus resulting in an Inexact Newton method. Here
right preconditioning is available by way of the preconditioning setup and solve
routines supplied by the user, in which case the iterative method is applied to
the linear systems :math:`(JP^{-1})(P\delta) = -F`, where :math:`P` denotes the
right preconditioning matrix.

Additionally, it is possible for users to supply a matrix-based iterative linear solver to KINSOL, resulting in a
Modified Inexact Newton method. As with the direct linear solvers, the Jacobian matrix is updated infrequently;
similarly as with iterative linear solvers the linear system is solved only approximately.

.. _KINSOL.Mathematics.JacUpdate:

Jacobian information update strategy
------------------------------------

In general, unless specified otherwise by the user, KINSOL strives to update Jacobian information (the actual
system Jacobian :math:`J` in the case of matrix-based linear solvers, and the preconditioner matrix :math:`P` in the
case of iterative linear solvers) as infrequently as possible to balance the high costs of matrix operations against
other costs. Specifically, these updates occur when:

-  the problem is initialized,

-  :math:`\|\lambda\delta_{n-1}\|_{D_u,\infty} > 1.5` (Inexact Newton only),

-  ``mbset``\ :math:`=10` nonlinear iterations have passed since the last update,

-  the linear solver failed recoverably with outdated Jacobian information,

-  the global strategy failed with outdated Jacobian information, or

-  :math:`\|\lambda\delta_{n}\|_{D_u,\infty} <` steptol with outdated Jacobian or preconditioner information,

where the norm :math:`\|\cdot\|_{D_u,\infty}` is defined below in :eq:`KIN_scaled-norms`.

KINSOL allows, through optional solver inputs, changes to the above strategy. Indeed, the user can disable the
initial Jacobian information evaluation or change the default value of ``mbset``, the number of nonlinear iterations
after which a Jacobian information update is enforced.

.. _KINSOL.Mathematics.Scaling:

Scaling
-------

To address the case of ill-conditioned nonlinear systems, KINSOL allows users to prescribe scaling factors both for the
solution vector and for the residual vector. For scaling to be used, the user should supply values :math:`D_u`, which
are diagonal elements of the scaling matrix such that :math:`D_u u_n` has all components roughly the same magnitude when
:math:`u_n` is close to a solution, and :math:`D_F`, which are diagonal scaling matrix elements such that :math:`D_F F`
has all components roughly the same magnitude when :math:`u_n` is not too close to a solution. Based on
these scaling matrices, we define the following scaled norms:

.. math::
   \|z\|_{D_u} = \|D_u z\|_2, \;\; \|z\|_{D_F} = \|D_F z\|_2, \;\;
   \|z\|_{D_u,\infty} = \|D_u z\|_\infty, \;\; {\rm and} \;\;
   \|z\|_{D_F,\infty} = \|D_F z\|_\infty
  :label: KIN_scaled-norms

where :math:`\| \cdot \|_\infty` is the max norm. When scaling values are provided for the solution
vector, these values are automatically incorporated into the calculation of the perturbations used for
the default difference quotient
approximations for Jacobian information; see :eq:`KIN_sigmaDQ_direct` and :eq:`KIN_sigmaDQ_iterative` below.

Globalization strategy
----------------------

Two methods of applying a computed step :math:`\delta_n` to the previously computed solution vector are implemented. The
first and simplest is the standard Newton strategy which applies step 2(b) as above with :math:`\lambda` always set to
:math:`1`. The other method is a global strategy, which attempts to use the direction implied by :math:`\delta_n` in the
most efficient way for furthering convergence of the nonlinear problem. This technique is implemented in the second
strategy, called Linesearch. This option employs both the :math:`\alpha` and :math:`\beta` conditions of the
Goldstein-Armijo linesearch given in :cite:p:`DeSc:96` for step 2(b), where :math:`\lambda` is chosen to
guarantee a sufficient decrease in :math:`F` relative to the step length as well as a minimum step length relative to
the initial rate of decrease of :math:`F`. One property of the algorithm is that the full Newton step tends to be taken
close to the solution.

KINSOL implements a backtracking algorithm to first find a value :math:`\lambda` such that
:math:`u_n + \lambda \delta_n` satisfies the sufficient decrease condition (or :math:`\alpha`-condition)

.. math:: F(u_n + \lambda\delta_n) \le F(u_n) + \alpha \nabla F(u_n)^T \lambda\delta_n \, ,

where :math:`\alpha = 10^{-4}`. Although backtracking in itself guarantees that the step is not too small, KINSOL
secondly relaxes :math:`\lambda` to satisfy the so-called :math:`\beta`-condition (equivalent to Wolfe’s curvature
condition):

.. math:: F(u_n + \lambda\delta_n) \ge F(u_n) + \beta \nabla F(u_n)^T \lambda\delta_n \, ,

where :math:`\beta = 0.9`. During this second phase, :math:`\lambda` is allowed to vary in the interval
:math:`[\lambda_{min} , \lambda_{max}]` where

.. math::
  \lambda_{min} =  \frac{{steptol}}{\| \bar\delta_n\|_\infty} \, , \quad
  \bar\delta_n^j = \frac{\delta_n^j}{1/D_u^j + |u^j|} \, ,

and :math:`\lambda_{max}` corresponds to the maximum feasible step size at the current iteration (typically
:math:`\lambda_{max} = {stepmax} / \|\delta_n\|_{D_u}`). In the above expressions, :math:`v^j` denotes the
:math:`j`\ th component of a vector :math:`v`.

For more details, the reader is referred to :cite:p:`DeSc:96`.

Nonlinear iteration stopping criteria
-------------------------------------

Stopping criteria for the Newton method are applied to both of the nonlinear residual and the step length. For the
former, the Newton iteration must pass a stopping test

.. math:: \|F(u_n)\|_{D_F,\infty} < \text{ftol} \, ,

where ftol is an input scalar tolerance with a default value of :math:`U^{1/3}`. Here :math:`U` is the machine unit
roundoff. For the latter, the Newton method will terminate when the maximum scaled step is below a given tolerance

.. math:: \|\lambda\delta_n\|_{D_u,\infty} < \text{steptol} \, ,

where steptol is an input scalar tolerance with a default value of :math:`U^{2/3}`. Only the first condition (small
residual) is considered a successful completion of KINSOL. The second condition (small step) may indicate that the
iteration is stalled near a point for which the residual is still unacceptable.

Additional constraints
----------------------

As a user option, KINSOL permits the application of inequality constraints, :math:`u^i > 0` and :math:`u^i < 0`,
as well as :math:`u^i \geq 0` and :math:`u^i \leq 0`, where :math:`u^i` is the :math:`i`\ th component of :math:`u`. Any
such constraint, or no constraint, may be imposed on each component. KINSOL will reduce step lengths in order to
ensure that no constraint is violated. Specifically, if a new Newton iterate will violate a constraint, the maximum step
length along the Newton direction that will satisfy all constraints is found, and :math:`\delta_n` in Step 2(b) is
scaled to take a step of that length.

.. _KINSOL.Mathematics.ModifiedNewtonResidualMon:

Residual monitoring for Modified Newton method
----------------------------------------------

When using a matrix-based linear solver, in addition to the strategy described above for the update of the Jacobian
matrix, KINSOL also provides an optional nonlinear residual monitoring scheme to control when the system Jacobian
is updated. Specifically, a Jacobian update will also occur when ``mbsetsub=5`` nonlinear iterations have
passed since the last update and

.. math:: \|F(u_n)\|_{D_F} > \omega \|F(u_m)\|_{D_F} \, ,

where :math:`u_n` is the current iterate and :math:`u_m` is the iterate at the last Jacobian update. The scalar
:math:`\omega` is given by

.. math::
   \omega = \min \left (\omega_{min} \, e^{\max \left ( 0 , \rho - 1 \right )} , \omega_{max}\right ) \, ,
   :label: KIN_resmon_omega

with :math:`\rho` defined as

.. math:: \rho = \frac{\|F(u_n) \|_{D_F}}{\text{ftol}} \, ,

where ftol is the input scalar tolerance discussed before. Optionally, a constant value :math:`\omega_{const}` can be
used for the parameter :math:`\omega`.

The constants controlling the nonlinear residual monitoring algorithm can be changed from their default values through
optional inputs to KINSOL. These include the parameters :math:`\omega_{min}` and :math:`\omega_{max}`, the
constant value :math:`\omega_{const}`, and the threshold ``mbsetsub``.

.. _KINSOL.Mathematics.InexactNewtonStopCrit:

Stopping criteria for iterative linear solvers
----------------------------------------------

When using an Inexact Newton method (i.e. when an iterative linear solver is used), the convergence of the overall
nonlinear solver is intimately coupled with the accuracy with which the linear solver in 2(a) above is solved.
KINSOL provides three options for stopping criteria for the linear system solver, including the two algorithms of
Eisenstat and Walker :cite:p:`EiWa:96`. More precisely, the Krylov iteration must pass a stopping test

.. math:: \|J \delta_n + F\|_{D_F} < (\eta_n + U) \|F\|_{D_F} \, ,

where :math:`\eta_n` is one of:

Eisenstat and Walker Choice 1
   .. math::

      \eta_n = \frac{\left|\; \|F(u_n)\|_{D_F}
            - \|F(u_{n-1}) + J(u_{n-1}) \delta_n \|_{D_F}
            \; \right|}
        {\|F(u_{n-1})\|_{D_F}} \, ,

Eisenstat and Walker Choice 2
   .. math::

      \eta_n = \gamma
        \left( \frac{ \|F(u_n)\|_{D_F}}{\|F(u_{n-1})\|_{D_F}} \right)^{\alpha} \, ,

   where default values of :math:`\gamma` and :math:`\alpha` are :math:`0.9` and :math:`2`, respectively.

Constant :math:`\eta`
   .. math:: \eta_n = \text{constant},

   with 0.1 as the default.

The default strategy is "Eisenstat and Walker Choice 1". For both options 1 and 2, appropriate safeguards are
incorporated to ensure that :math:`\eta` does not decrease too quickly :cite:p:`EiWa:96`.

Difference quotient Jacobian approximations
-------------------------------------------

With the :ref:`SUNMATRIX_DENSE <SUNMatrix.Dense>` and :ref:`SUNMATRIX_BAND <SUNMatrix.Band>` matrix modules,
the Jacobian may be supplied by a user routine, or approximated
by difference quotients, at the user’s option. In the latter case, we use the usual approximation

.. math::
   J^{ij} = [F^i(u+\sigma_j e^j) - F^i(u)]/\sigma_j \, .
   :label: KIN_JacDQ

The increments :math:`\sigma_j` are given by

.. math::
   \sigma_j = \sqrt{U} \; \max\left\{ |u^j| , 1/D_u^j \right\} \, .
   :label: KIN_sigmaDQ_direct

In the dense case, this scheme requires :math:`N` evaluations of :math:`F`, one for each column of :math:`J`. In the
band case, the columns of :math:`J` are computed in groups, by the Curtis-Powell-Reid algorithm, with the number of
:math:`F` evaluations equal to the bandwidth. The parameter :math:`U` above can (optionally) be replaced by a
user-specified value, ``relfunc``.

We note that with sparse and user-supplied matrix-based linear solvers, the Jacobian *must* be supplied by a user
routine, i.e. it is not approximated internally within KINSOL.

In the case of a matrix-free iterative linear solver, Jacobian information is needed only as matrix-vector products
:math:`Jv`. If a routine for :math:`Jv` is not supplied, these products are approximated by directional difference
quotients as

.. math::
   J(u) v \approx [F(u+\sigma v) - F(u)]/\sigma \, ,
   :label: KIN_JvDQ

where :math:`u` is the current approximation to a root of :eq:`KIN_nonlinear_system`, and
:math:`\sigma` is a scalar. The choice of :math:`\sigma` is taken from :cite:p:`BrSa:90` and is given by

.. math::
   \sigma = \frac{\max \{|u^T v|, u^T_{typ} |v|\}}{\|v\|_2^2}
   \mbox{sign}(u^T v) \sqrt{U} \, ,
   :label: KIN_sigmaDQ_iterative

where :math:`u_{typ}` is a vector of typical values for the absolute values of the solution (and can be taken to be
inverses of the scale factors given for :math:`u` as described below). This formula is suitable for *scaled* vectors
:math:`u` and :math:`v`, and so is applied to :math:`D_u u` and :math:`D_u v`. The parameter :math:`U` above can
(optionally) be replaced by a user-specified value, ``relfunc``. Convergence of the Newton method is maintained as long
as the value of :math:`\sigma` remains appropriately small, as shown in :cite:p:`Bro:87`.

Basic Fixed Point iteration
---------------------------

The basic fixed-point iteration scheme implemented in KINSOL is given by:

#. Set :math:`u_0 =` an initial guess

#. For :math:`n = 0, 1, 2,...` until convergence do:

   -  Set :math:`u_{n+1} = (1 - \beta) u_n + \beta G(u_n)`.

   -  Test for convergence.

Here, :math:`u_n` is the :math:`n`-th iterate to :math:`u`. At each stage in the iteration process, the function
:math:`G` is applied to the current iterate with the damping parameter :math:`\beta` to produce a new iterate,
:math:`u_{n+1}`. A test for convergence is made before the iteration continues.

For Picard iteration, as implemented in KINSOL, we consider a special form of the nonlinear function :math:`F`,
such that :math:`F(u) = Lu - N(u)`, where :math:`L` is a constant nonsingular matrix and :math:`N` is (in general)
nonlinear. Then the fixed-point function :math:`G` is defined as :math:`G(u) = u - L^{-1}F(u)`. The Picard iteration is
given by:

#. Set :math:`u_0 =` an initial guess

#. For :math:`n = 0, 1, 2,...` until convergence do:

   -  Set :math:`u_{n+1} = (1 - \beta) u_n + \beta G(u_n)` where :math:`G(u_n) \equiv u_n - L^{-1}F(u_n)`.

   -  Test :math:`F(u_{n+1})` for convergence.

Here, :math:`u_n` is the :math:`n`-th iterate to :math:`u`. Within each iteration, the Picard step is computed then
added to :math:`u_n` with the damping parameter :math:`\beta` to produce the new iterate. Next, the nonlinear residual
function is evaluated at the new iterate, and convergence is checked. Noting that :math:`L^{-1}N(u) = u - L^{-1}F(u)`,
the above iteration can be written in the same form as a Newton iteration except that here, :math:`L` is in the role of
the Jacobian. Within KINSOL, however, we leave this in a fixed-point form as above. For more information,
see page 182 of :cite:p:`Ortega-Rheinbolt00`.

Anderson Acceleration
---------------------

The Picard and fixed point methods can be significantly accelerated using Anderson’s method
:cite:p:`Anderson65, Walker-Ni09, Fang-Saad09, LWWY11`. Anderson acceleration can be formulated as follows:

1. Set :math:`u_0 =` an initial guess and :math:`m \ge 1`

2. Set :math:`u_1 = G(u_0)`

3. For :math:`n = 1, 2,...` until convergence do:

   a.  Set :math:`m_n = \min\{m,n\}`

   b.  Set :math:`F_{n} = (f_{n-m_n}, \ldots, f_n)`, where :math:`f_i=G(u_i)-u_i`

   c.  Determine :math:`\alpha^{(n)} = (\alpha_0^{(n)}, \ldots, \alpha_{m_n}^{(n)})` that solves
       :math:`\displaystyle\min_\alpha  \| F_n \alpha^T \|_2` such that :math:`\displaystyle\sum_{i=0}^{m_n} \alpha_i = 1`

   d.  Set :math:`\displaystyle u_{n+1} = \beta \sum_{i=0}^{m_n} \alpha_i^{(n)} G(u_{n-m_n+i}) + (1-\beta) \sum_{i=0}^{m_n} \alpha_i^{(n)} u_{n-m_{n}+i}`

   e.  Test for convergence

It has been implemented in KINSOL by turning the constrained linear least-squares problem in step 3c into an
unconstrained one leading to the algorithm given below:

1. Set :math:`u_0 =` an initial guess and :math:`m \ge 1`

2. Set :math:`u_1 = G(u_0)`

3. For :math:`n = 1, 2,...` until convergence do:

   a. Set :math:`m_n = \min\{m,n\}`

   b. Set :math:`\Delta F_{n} = (\Delta f_{n-m_n}, \ldots, \Delta f_{n-1})`, where :math:`\Delta f_i = f_{i+1} - f_i`
      and :math:`f_i=G(u_i)-u_i`

   c. Determine :math:`\gamma^{(n)} = (\gamma_0^{(n)}, \ldots, \gamma_{m_n-1}^{(n)})` that solves
      :math:`\displaystyle\min_\gamma  \| f_n - \Delta F_n \gamma^T \|_2`

   d. Set :math:`\displaystyle u_{n+1} = G(u_n)-\sum_{i=0}^{m_n-1} \gamma_i^{(n)} \Delta g_{n-m_n+i} - (1-\beta)(f(u_n) - \sum_{i=0}^{m_n-1} \gamma_i^{(n)} \Delta f_{n-m_n+i})` with
      :math:`\Delta g_i = G(u_{i+1}) - G(u_i)`

   e. Test for convergence

The least-squares problem in 3c is solved by applying a QR factorization to :math:`\Delta F_n = Q_n R_n` and solving
:math:`R_n \gamma = Q_n^T f_n`. By default the damping is disabled i.e., :math:`\beta = 1.0`.

The Anderson acceleration implementation includes an option to delay the start of acceleration until after a given
number of initial fixed-point or Picard iterations have been completed. This delay can be beneficial when the underlying
method has strong global convergence properties as the initial iterations may help bring the iterates closer to a
solution before starting the acceleration.

.. _Anderson_QR:

Anderson Acceleration QR Factorization
--------------------------------------

The default QR factorization routine used in Anderson acceleration is Modified
Gram-Schmidt, a stable orthogonalization routine that requires an increasing
number of synchronizations per iteration dependent upon the number of vectors
being orthgonalized against. While practical use of Anderson acceleration only
requires a small number of vectors to be used in the QR factorization, this
linearly scaling number of synchronizations per iteration can yield poor
performance when Anderson acceleration is performed in a parallel setting. To
combat this poor performance, low synchronization QR routines are available to
the user, in particular: Inverse Compact WY Modified Gram-Schmidt
:cite:p:`lowSyncGMRES`, along with variants of Classical Gram-Schmidt with
Reorthogonalization :cite:p:`hernandez2005parallel`.  While all of these QR
factorization routines are mathematically equivalent, they do not exhibit the
same stability when performed with floating point arithmetic or in a parallel
setting.

Inverse Compact WY Modified Gram-Schmidt, which is based on triangular solve
variants of Gram-Schmidt that were developed within the context of GMRES, is an
option that only requires two synchronizations per iteration. Additionally, it
adds a lower triangular solve at every iteration, but this generally does not
affect performance due to the system solve being small i.e., the number of
vectors being orthgonalized against.

The remaining orthogonalization options are based on and include Classical
Gram-Schmidt with Reorthogonalization (CGS-2). CGS-2 only requires three
synchronizations per iteration, but does not exhibit the same stability as
Modified Gram-Schmidt. Classical Gram-Schmidt with Delayed Reorthogolonization
has the same stability as CGS-2, but it reduces the number of synchronizations
per iteration to two.

Fixed-point - Anderson Acceleration Stopping Criterion
------------------------------------------------------

The default stopping criterion is

.. math:: \|u_{n+1} - u_{n} \|_{D_F,\infty} < \text{gtol} \, ,

where :math:`D_F` is a user-defined diagonal matrix that can be the identity or a scaling matrix chosen so that the
components of :math:`D_F (G(u)-u)` have roughly the same order of magnitude. Note that when using Anderson acceleration,
convergence is checked after the acceleration is applied.

Picard - Anderson Acceleration Stopping Criterion
-------------------------------------------------

The default stopping criterion is

.. math:: \|F(u_{n+1})\|_{D_F,\infty} < \text{ftol} \, ,

where :math:`D_F` is a user-defined diagonal matrix that can be the identity or a scaling matrix chosen so that the
components of :math:`D_F F(u)` have roughly the same order of magnitude. Note that when using Anderson acceleration,
convergence is checked after the acceleration is applied.
