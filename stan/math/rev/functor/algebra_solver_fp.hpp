#ifndef STAN_MATH_REV_FUNCTOR_FP_SOLVER_HPP
#define STAN_MATH_REV_FUNCTOR_FP_SOLVER_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/fun/value_of.hpp>
#include <stan/math/rev/functor/algebra_system.hpp>
#include <stan/math/rev/functor/kinsol_data.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/mdivide_left.hpp>
#include <stan/math/prim/fun/to_array_1d.hpp>
#include <stan/math/prim/fun/to_vector.hpp>
#include <stan/math/prim/fun/value_of.hpp>
#include <sundials/sundials_context.h>
#include <kinsol/kinsol.h>
#include <sunmatrix/sunmatrix_dense.h>
#include <sunlinsol/sunlinsol_dense.h>
#include <nvector/nvector_serial.h>
#include <algorithm>
#include <iostream>
#include <string>
#include <vector>

namespace stan {
namespace math {

/**
 * KINSOL algebraic system data holder that handles
 * construction & destruction of SUNDIALS data, as well as
 * auxiliary data that will be used for functor evaluation.
 *
 * @tparam F functor type for system function.
 */
template <typename F>
struct KinsolFixedPointEnv {
  /** Sundials context **/
  sundials::Context sundials_context_;
  /** RHS functor. */
  const F& f_;
  /** val of params for @c y_ to refer to when
   params are @c var type */
  const Eigen::VectorXd y_dummy;
  /** ref to val of params */
  const Eigen::VectorXd& y_;
  /** system size */
  const size_t N_;
  /** nb. of params */
  const size_t M_;
  /** real data */
  const std::vector<double>& x_r_;
  /** integer data */
  const std::vector<int>& x_i_;
  /** message stream */
  std::ostream* msgs_;
  /** KINSOL memory block */
  void* mem_;
  /** NVECTOR for unknowns */
  N_Vector nv_x_;
  /** NVECTOR for scaling u */
  N_Vector nv_u_scal_;
  /** NVECTOR for scaling f */
  N_Vector nv_f_scal_;

  /** Constructor when y is data */
  template <typename T, typename T_u, typename T_f>
  KinsolFixedPointEnv(const F& f, const Eigen::Matrix<T, -1, 1>& x,
                      const Eigen::VectorXd& y, const std::vector<double>& x_r,
                      const std::vector<int>& x_i, std::ostream* msgs,
                      const std::vector<T_u>& u_scale,
                      const std::vector<T_f>& f_scale)
      : sundials_context_(),
        f_(f),
        y_dummy(),
        y_(y),
        N_(x.size()),
        M_(y.size()),
        x_r_(x_r),
        x_i_(x_i),
        msgs_(msgs),
        mem_(KINCreate(sundials_context_)),
        nv_x_(N_VNew_Serial(N_, sundials_context_)),
        nv_u_scal_(N_VNew_Serial(N_, sundials_context_)),
        nv_f_scal_(N_VNew_Serial(N_, sundials_context_)) {
    for (int i = 0; i < N_; ++i) {
      NV_Ith_S(nv_x_, i) = stan::math::value_of(x(i));
      NV_Ith_S(nv_u_scal_, i) = stan::math::value_of(u_scale[i]);
      NV_Ith_S(nv_f_scal_, i) = stan::math::value_of(f_scale[i]);
    }
  }

  /** Constructor when y is param */
  template <typename T, typename T_u, typename T_f>
  KinsolFixedPointEnv(const F& f, const Eigen::Matrix<T, -1, 1>& x,
                      const Eigen::Matrix<stan::math::var, -1, 1>& y,
                      const std::vector<double>& x_r,
                      const std::vector<int>& x_i, std::ostream* msgs,
                      const std::vector<T_u>& u_scale,
                      const std::vector<T_f>& f_scale)
      : sundials_context_(),
        f_(f),
        y_dummy(stan::math::value_of(y)),
        y_(y_dummy),
        N_(x.size()),
        M_(y.size()),
        x_r_(x_r),
        x_i_(x_i),
        msgs_(msgs),
        mem_(KINCreate(sundials_context_)),
        nv_x_(N_VNew_Serial(N_, sundials_context_)),
        nv_u_scal_(N_VNew_Serial(N_, sundials_context_)),
        nv_f_scal_(N_VNew_Serial(N_, sundials_context_)) {
    for (int i = 0; i < N_; ++i) {
      NV_Ith_S(nv_x_, i) = stan::math::value_of(x(i));
      NV_Ith_S(nv_u_scal_, i) = stan::math::value_of(u_scale[i]);
      NV_Ith_S(nv_f_scal_, i) = stan::math::value_of(f_scale[i]);
    }
  }

  ~KinsolFixedPointEnv() {
    N_VDestroy_Serial(nv_x_);
    N_VDestroy_Serial(nv_u_scal_);
    N_VDestroy_Serial(nv_f_scal_);
    KINFree(&mem_);
  }

  /** Implements the user-defined function passed to KINSOL. */
  static int kinsol_f_system(N_Vector x, N_Vector f, void* user_data) {
    auto g = static_cast<const KinsolFixedPointEnv<F>*>(user_data);
    Eigen::VectorXd x_eigen(Eigen::Map<Eigen::VectorXd>(NV_DATA_S(x), g->N_));
    Eigen::Map<Eigen::VectorXd>(N_VGetArrayPointer(f), g->N_)
        = g->f_(x_eigen, g->y_, g->x_r_, g->x_i_, g->msgs_);
    return 0;
  }
};

/**
 * Calculate Jacobian Jxy(Jacobian of unknown x w.r.t. the * param y)
 * given the solution. Specifically, for
 *
 *  x - f(x, y) = 0
 *
 * we have (Jpq = Jacobian matrix dq/dq)
 *
 * Jxy - Jfx * Jxy = Jfy
 *
 * therefore Jxy can be solved from system
 *
 * (I - Jfx) * Jxy = Jfy
 *
 * Jfx and Jfy are obtained through one AD evaluation of f
 * w.r.t combined vector [x, y].
 */
struct FixedPointADJac {
  /**
   * Calculate Jacobian Jxy.
   *
   * @tparam F RHS functor type
   * @param x fixed point solution
   * @param y RHS parameters
   * @param env KINSOL working environment, see doc for @c KinsolFixedPointEnv.
   */
  template <typename F>
  inline Eigen::Matrix<stan::math::var, -1, 1> operator()(
      const Eigen::VectorXd& x, const Eigen::Matrix<stan::math::var, -1, 1>& y,
      KinsolFixedPointEnv<F>& env) {
    using stan::math::precomputed_gradients;
    using stan::math::to_array_1d;
    using stan::math::var;

    auto g = [&env](const Eigen::Matrix<var, -1, 1>& par_) {
      Eigen::Matrix<var, -1, 1> x_(par_.head(env.N_));
      Eigen::Matrix<var, -1, 1> y_(par_.tail(env.M_));
      return env.f_(x_, y_, env.x_r_, env.x_i_, env.msgs_);
    };

    Eigen::VectorXd theta(x.size() + y.size());
    for (int i = 0; i < env.N_; ++i) {
      theta(i) = x(i);
    }
    for (int i = 0; i < env.M_; ++i) {
      theta(i + env.N_) = env.y_(i);
    }
    Eigen::Matrix<double, -1, 1> fx;
    Eigen::Matrix<double, -1, -1> J_theta;
    stan::math::jacobian(g, theta, fx, J_theta);
    Eigen::MatrixXd A(J_theta.block(0, 0, env.N_, env.N_));
    Eigen::MatrixXd b(J_theta.block(0, env.N_, env.N_, env.M_));
    A = Eigen::MatrixXd::Identity(env.N_, env.N_) - A;
    Eigen::MatrixXd Jxy = A.colPivHouseholderQr().solve(b);
    std::vector<double> gradients(env.M_);
    Eigen::Matrix<var, -1, 1> x_sol(env.N_);
    std::vector<stan::math::var> yv(to_array_1d(y));
    for (int i = 0; i < env.N_; ++i) {
      gradients = to_array_1d(Eigen::VectorXd(Jxy.row(i)));
      x_sol[i] = precomputed_gradients(x(i), yv, gradients);
    }
    return x_sol;
  }
};

/**
 * Fixed point solver for problem of form
 *
 * x = F(x; theta)
 *
 * with x as unknowns and theta parameters.
 *
 * The solution for FP iteration
 * doesn't involve derivatives but only data types.
 *
 * @tparam fp_env_type solver environment setup that handles
 *                     workspace & auxiliary data encapsulation & RAII, namely
 *                     the work environment. Currently only support KINSOL's
 *                     dense matrix.
 * @tparam fp_jac_type functor type for calculating the
 *                     jacobian. Currently only support @c
 *                     FixedPointADJac that obtain dense Jacobian
 *                     through QR decomposition.
 */
template <typename fp_env_type, typename fp_jac_type>
struct FixedPointSolver;

/**
 * Specialization for fixed point solver when using KINSOL.
 *
 * @tparam F RHS functor for fixed point iteration.
 * @tparam fp_jac_type functor type for calculating the jacobian
 */
template <typename F, typename fp_jac_type>
struct FixedPointSolver<KinsolFixedPointEnv<F>, fp_jac_type> {
  /**
   * Solve FP using KINSOL
   *
   * @param x initial point and final solution.
   * @param env KINSOL solution environment
   * @param f_tol Function tolerance
   * @param max_num_steps max nb. of iterations.
   */
  void kinsol_solve_fp(Eigen::VectorXd& x, KinsolFixedPointEnv<F>& env,
                       double f_tol, int max_num_steps) {
    int N = env.N_;
    void* mem = env.mem_;

    const int default_anderson_depth = 4;
    int anderson_depth = std::min(N, default_anderson_depth);

    CHECK_KINSOL_CALL(KINSetNumMaxIters(mem, max_num_steps));
    CHECK_KINSOL_CALL(KINSetMAA(mem, anderson_depth));
    CHECK_KINSOL_CALL(KINInit(mem, &env.kinsol_f_system, env.nv_x_));
    CHECK_KINSOL_CALL(KINSetFuncNormTol(mem, f_tol));
    CHECK_KINSOL_CALL(KINSetUserData(mem, static_cast<void*>(&env)));
    kinsol_check(KINSol(mem, env.nv_x_, KIN_FP, env.nv_u_scal_, env.nv_f_scal_),
                 "KINSol", max_num_steps);

    for (int i = 0; i < N; ++i) {
      x(i) = NV_Ith_S(env.nv_x_, i);
    }
  }

  /**
   * Solve data-only FP problem so no need to calculate jacobian.
   *
   * @tparam T1 type of init point of iterations
   *
   * @param x initial point and final solution.
   * @param y RHS functor parameters
   * @param env KINSOL solution environment
   * @param f_tol Function tolerance
   * @param max_num_steps max nb. of iterations.
   */
  template <typename T1>
  Eigen::Matrix<double, -1, 1> solve(const Eigen::Matrix<T1, -1, 1>& x,
                                     const Eigen::Matrix<double, -1, 1>& y,
                                     KinsolFixedPointEnv<F>& env, double f_tol,
                                     int max_num_steps) {
    Eigen::VectorXd xd(stan::math::value_of(x));
    kinsol_solve_fp(xd, env, f_tol, max_num_steps);
    return xd;
  }

  /**
   * Solve FP problem and calculate jacobian.
   *
   * @tparam T1 type of init point of iterations
   *
   * @param x initial point and final solution.
   * @param y RHS functor parameters
   * @param env KINSOL solution environment
   * @param f_tol Function tolerance
   * @param max_num_steps max nb. of iterations.
   */
  template <typename T1>
  Eigen::Matrix<stan::math::var, -1, 1> solve(
      const Eigen::Matrix<T1, -1, 1>& x,
      const Eigen::Matrix<stan::math::var, -1, 1>& y,
      KinsolFixedPointEnv<F>& env, double f_tol, int max_num_steps) {
    using stan::math::value_of;
    using stan::math::var;

    // FP solution
    Eigen::VectorXd xd(solve(x, Eigen::VectorXd(), env, f_tol, max_num_steps));

    fp_jac_type jac_sol;
    return jac_sol(xd, y, env);
  }
};

/**
 * Return a fixed pointer to the specified system of algebraic
 * equations of form
 *
 * x = F(x; theta)
 *
 * given an initial guess x, and parameters theta and data. Use the
 * KINSOL solver from the SUNDIALS suite.
 *
 * The user can also specify the scaling controls, the function
 * tolerance, and the maximum number of steps.
 *
 * @tparam F type of equation system function.
 * @tparam T type of initial guess vector. The final solution
 *           type doesn't depend on initial guess type,
 *           but we allow initial guess to be either data or param.
 * @tparam T_u type of scaling vector for unknowns. We allow
 *             it to be @c var because scaling could be parameter
 *             dependent. Internally these params are converted to data
 *             because scaling is applied.
 * @tparam T_f type of scaling vector for residual. We allow
 *             it to be @c var because scaling could be parameter
 *             dependent. Internally these params are converted to data
 *             because scaling is applied.
 *
 * @param[in] f Functor that evaluated the system of equations.
 * @param[in] x Vector of starting values.
 * @param[in] y Parameter vector for the equation system. The function
 *            is overloaded to treat y as a vector of doubles or of a
 *            a template type T.
 * @param[in] x_r Continuous data vector for the equation system.
 * @param[in] x_i Integer data vector for the equation system.
 * @param[in, out] msgs The print stream for warning messages.
 * @param[in] u_scale diagonal scaling matrix elements Du
 *                    such that Du*x has all components roughly the same
 *                    magnitude when x is close to a solution.
 *                    (ref. KINSOL user guide chap.2 sec. "Scaling")
 * @param[in] f_scale diagonal scaling matrix elements such
 *                    that Df*(x-f(x)) has all components roughly the same
 *                    magnitude when x is not too close to a solution.
 *                    (ref. KINSOL user guide chap.2 sec. "Scaling")
 * @param[in] f_tol Function-norm stopping tolerance.
 * @param[in] max_num_steps maximum number of function evaluations.
 * @throw <code>std::invalid_argument</code> if x has size zero.
 * @throw <code>std::invalid_argument</code> if x has non-finite elements.
 * @throw <code>std::invalid_argument</code> if y has non-finite elements.
 * @throw <code>std::invalid_argument</code> if dat has non-finite elements.
 * @throw <code>std::invalid_argument</code> if dat_int has non-finite elements.
 * @throw <code>std::invalid_argument</code> if scaled_step_size is strictly
 * negative.
 * @throw <code>std::invalid_argument</code> if function_tolerance is strictly
 * negative.
 * @throw <code>std::invalid_argument</code> if max_num_steps is not positive.
 * @throw <code>boost::math::evaluation_error</code> (which is a subclass of
 * <code>std::runtime_error</code>) if solver exceeds max_num_steps.
 */
template <typename F, typename T1, typename T2, typename T_u, typename T_f>
Eigen::Matrix<T2, -1, 1> algebra_solver_fp(
    const F& f, const Eigen::Matrix<T1, -1, 1>& x,
    const Eigen::Matrix<T2, -1, 1>& y, const std::vector<double>& x_r,
    const std::vector<int>& x_i, const std::vector<T_u>& u_scale,
    const std::vector<T_f>& f_scale, std::ostream* msgs = nullptr,
    double f_tol = 1e-8,
    int max_num_steps = 200) {  // NOLINT(runtime/int)
  algebra_solver_check(x, y, x_r, x_i, f_tol, max_num_steps);
  check_nonnegative("algebra_solver", "u_scale", u_scale);
  check_nonnegative("algebra_solver", "f_scale", f_scale);
  check_matching_sizes("algebra_solver", "the algebraic system's output",
                       value_of(f(x, y, x_r, x_i, msgs)),
                       "the vector of unknowns, x,", x);

  KinsolFixedPointEnv<F> env(f, x, y, x_r, x_i, msgs, u_scale,
                             f_scale);  // NOLINT
  FixedPointSolver<KinsolFixedPointEnv<F>, FixedPointADJac> fp;
  return fp.solve(x, y, env, f_tol, max_num_steps);
}

}  // namespace math
}  // namespace stan

#endif
