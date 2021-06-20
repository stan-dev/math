#ifndef STAN_MATH_REV_FUNCTOR_FP_SOLVER_HPP
#define STAN_MATH_REV_FUNCTOR_FP_SOLVER_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/fun/value_of.hpp>
#include <stan/math/rev/functor/algebra_system.hpp>
#include <stan/math/rev/functor/kinsol_data.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/value_of.hpp>
#include <stan/math/prim/fun/eval.hpp>
#include <stan/math/prim/functor/apply.hpp>
#include <stan/math/prim/functor/algebra_solver_adapter.hpp>
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
template <typename F, typename T_u, typename T_f, typename... Args>
struct KinsolFixedPointEnv {
  /** RHS functor. */
  const F& f_;
  /** system size */
  const size_t N_;
  /** message stream */
  std::ostream* msgs_;
  /** arguments and parameters */
  std::tuple<const Args&...> args_tuple_;
  /** KINSOL memory block */
  void* mem_;
  /** NVECTOR for unknowns */
  N_Vector nv_x_;
  /** NVECTOR for scaling u */
  N_Vector nv_u_scal_;
  /** NVECTOR for scaling f */
  N_Vector nv_f_scal_;

  KinsolFixedPointEnv(const F& f, const Eigen::MatrixXd x,
                      const std::vector<T_u>& u_scale,
                      const std::vector<T_f>& f_scale,
                      std::ostream* msgs, const Args&... args)
      : f_(f),
        N_(x.size()),
        msgs_(msgs),
        args_tuple_(args...),
        mem_(KINCreate()),
        nv_x_(N_VNew_Serial(N_)),
        nv_u_scal_(N_VNew_Serial(N_)),
        nv_f_scal_(N_VNew_Serial(N_)) {
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
    auto g =
      static_cast<const KinsolFixedPointEnv<F, T_u, T_f, Args...>*>(user_data);
    Eigen::VectorXd x_eigen(Eigen::Map<Eigen::VectorXd>(NV_DATA_S(x), g->N_));
    Eigen::Map<Eigen::VectorXd> f_map(N_VGetArrayPointer(f), g->N_);

    f_map = apply(
        [&](const auto&... args) {
          return g->f_(x_eigen, g->msgs_, args...);
        },
        g->args_tuple_);
    return 0;
  }
};

/**
 * Solve FP using KINSOL
 *
 * @param x initial point and final solution.
 * @param env KINSOL solution environment
 * @param f_tol Function tolerance
 * @param max_num_steps max nb. of iterations.
 */
template <typename F, typename T_u, typename T_f, typename... Args>
Eigen::VectorXd kinsol_solve_fp(
    const F& f,
    const Eigen::VectorXd& x,
    double function_tolerance,
    double max_num_steps,
    const std::vector<T_u>& u_scale,
    const std::vector<T_f>& f_scale,
    std::ostream* msgs,
    const Args&... args) {
  KinsolFixedPointEnv<F, T_u, T_f, Args...> env(
      f, x, u_scale, f_scale, msgs, args...);
  int N = env.N_;
  void* mem = env.mem_;
  const int default_anderson_depth = 4;
  int anderson_depth = std::min(N, default_anderson_depth);
  Eigen::VectorXd x_solution(N);

  check_flag_sundials(KINSetNumMaxIters(mem, max_num_steps),
                      "KINSetNumMaxIters");
  check_flag_sundials(KINSetMAA(mem, anderson_depth), "KINSetMAA");
  check_flag_sundials(KINInit(mem, &env.kinsol_f_system, env.nv_x_),
                      "KINInit");
  check_flag_sundials(KINSetFuncNormTol(mem, function_tolerance),
                      "KINSetFuncNormTol");
  check_flag_sundials(KINSetUserData(mem, static_cast<void*>(&env)),
                      "KINSetUserData");

  check_flag_kinsol(
      KINSol(mem, env.nv_x_, KIN_FP, env.nv_u_scal_, env.nv_f_scal_),
      max_num_steps);

  for (int i = 0; i < N; i++)
    x_solution(i) = NV_Ith_S(env.nv_x_, i);

  return x_solution;
}

/** Implementation of ordinary fixed point solver. */
template <typename F, typename T, typename T_u, typename T_f,
          typename... Args,
          require_eigen_vector_t<T>* = nullptr,
          require_all_st_arithmetic<Args...>* = nullptr>
Eigen::VectorXd algebra_solver_fp_impl(
    const F& f, const T& x, const double function_tolerance,
    const int max_num_steps, const std::vector<T_u>& u_scale,
    const std::vector<T_f>& f_scale, std::ostream* msgs, const Args&... args) {
  const auto& x_ref = to_ref(value_of(x));

  check_nonzero_size("algebra_solver_fp", "initial guess", x_ref);
  check_finite("algebra_solver_fp", "initial guess", x_ref);
  check_nonnegative("algebra_solver_fp", "u_scale", u_scale);
  check_nonnegative("algebra_solver_fp", "f_scale", f_scale);
  check_nonnegative("algebra_solver_fp", "function_tolerance",
                    function_tolerance);
  check_positive("algebra_solver_fp", "max_num_steps", max_num_steps);
  check_matching_sizes("algebra_solver_fp", "the algebraic system's output",
                       value_of(f(x_ref, msgs, args...)),
                       "the vector of unknowns, x,", x_ref);

  return kinsol_solve_fp(f, x_ref, function_tolerance, max_num_steps, u_scale,
                         f_scale, msgs, args...);
}

/* Implementation of autodiff fixed point solver. The Jacobian Jxy(Jacobian of
 * unknown x w.r.t. the * param y) is calculated given the solution as follows.
 * Since
 *
 *  x - f(x, y) = 0
 *
 * we have (Jpq being the Jacobian matrix dq/dq)
 *
 * Jxy - Jfx * Jxy = Jfy
 *
 * therefore Jxy can be solved from system
 *
 * (I - Jfx) * Jxy = Jfy
 *
 * Let eta be the adjoint with respect to x; then to calculate
 *
 * eta * Jxy
 *
 * we solve
 *
 * (eta * (I - Jfx)^(-1)) * Jfy
 *
 * (This is virtually identical to the Powell and Newton solvers, except Jfx
 * has been replaced by I - Jfx.)
 */
template <typename F, typename T, typename T_u, typename T_f,
          typename... Args,
          require_eigen_vector_t<T>* = nullptr,
          require_any_st_var<Args...>* = nullptr>
Eigen::Matrix<var, Eigen::Dynamic, 1> algebra_solver_fp_impl(
    const F& f, const T& x, const double function_tolerance,
    const int max_num_steps, const std::vector<T_u>& u_scale,
    const std::vector<T_f>& f_scale, std::ostream* msgs, const Args&... args) {
  const auto& x_val = to_ref(value_of(x));
  auto arena_args_tuple = std::make_tuple(to_arena(args)...);
  auto args_vals_tuple = std::make_tuple(eval(value_of(args))...);

  auto f_wrt_x = [&](const auto& x) {
    return apply([&](const auto&... args) { return f(x, msgs, args...); },
                 args_vals_tuple);
  };

  // FP solution
  Eigen::VectorXd theta_dbl = apply(
    [&](const auto&... vals) {
      return kinsol_solve_fp(f, x_val, function_tolerance, max_num_steps,
                             u_scale, f_scale, msgs, vals...);
    },
    args_vals_tuple);

  Eigen::MatrixXd Jf_x;
  Eigen::VectorXd f_x;

  jacobian(f_wrt_x, theta_dbl, f_x, Jf_x);
  int N = x.size();
  Jf_x = Eigen::MatrixXd::Identity(x.size(), x.size()) - Jf_x;

  using ret_type = Eigen::Matrix<var, Eigen::Dynamic, -1>;
  auto arena_Jf_x = to_arena(Jf_x);

  arena_t<ret_type> ret = theta_dbl;

  reverse_pass_callback([f, ret, arena_args_tuple, arena_Jf_x, msgs]() mutable {
    using Eigen::Dynamic;
    using Eigen::Matrix;
    using Eigen::MatrixXd;
    using Eigen::VectorXd;

    // Contract specificities with inverse Jacobian of f with respect to x.
    VectorXd ret_adj = ret.adj();
    VectorXd eta = arena_Jf_x.transpose().lu().solve(ret_adj);

    // Contract with Jacobian of f with respect to y using a nested reverse
    // autodiff pass.
    {
      nested_rev_autodiff rev;

      VectorXd ret_val = ret.val();
      auto x_nrad_ = apply(
          [&](const auto&... args) { return eval(f(ret_val, msgs, args...)); },
          arena_args_tuple);
      x_nrad_.adj() = eta;
      grad();
    }
  });

  return ret_type(ret);
}


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
    const Eigen::Matrix<T2, -1, 1>& y, const std::vector<double>& dat,
    const std::vector<int>& dat_int, const std::vector<T_u>& u_scale,
    const std::vector<T_f>& f_scale, std::ostream* msgs = nullptr,
    double function_tolerance = 1e-8,
    int max_num_steps = 200) {  // NOLINT(runtime/int)
  return algebra_solver_fp_impl(algebra_solver_adapter<F>(f), x, function_tolerance,
                           max_num_steps, u_scale, f_scale, msgs, y, dat,
                           dat_int);
}

}  // namespace math
}  // namespace stan

#endif
