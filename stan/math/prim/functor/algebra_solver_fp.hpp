#ifndef STAN_MATH_PRIM_FUNCTOR_FP_SOLVER_HPP
#define STAN_MATH_PRIM_FUNCTOR_FP_SOLVER_HPP

#include <stan/math/prim/meta.hpp>
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
 * @tparam T_u type of scaling vector for unknowns. We allow
 *             it to be @c var because scaling could be parameter
 *             dependent. Internally these params are converted to data
 *             because scaling is applied.
 * @tparam T_f type of scaling vector for residual. We allow
 *             it to be @c var because scaling could be parameter
 *             dependent. Internally these params are converted to data
 *             because scaling is applied.
 * @tparam Args types of additional parameters to the equation system functor.
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
                      const std::vector<T_f>& f_scale, std::ostream* const msgs,
                      const Args&... args)
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
    }
    for (int i = 0; i < N_; ++i) {
      NV_Ith_S(nv_u_scal_, i) = stan::math::value_of(u_scale[i]);
    }
    for (int i = 0; i < N_; ++i) {
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
  static int kinsol_f_system(const N_Vector x, const N_Vector f,
                             void* const user_data) {
    auto g = static_cast<const KinsolFixedPointEnv<F, T_u, T_f, Args...>*>(
        user_data);
    Eigen::Map<Eigen::VectorXd>(N_VGetArrayPointer(f), g->N_) = apply(
        [&x, &g](const auto&... args) {
          return g->f_(Eigen::Map<Eigen::VectorXd>(NV_DATA_S(x), g->N_), g->msgs_, args...);
        },
        g->args_tuple_);
    return 0;
  }
};

/**
 * Private interface for solving fixed point problems using KINSOL. Users should
 * call the KINSOL fixed point solver through `algebra_solver_fp` or
 * `algebra_solver_fp_impl`.
 *
 * @tparam F type of the equation system functor f
 * @tparam T_u type of scaling vector for unknowns. We allow
 *             it to be @c var because scaling could be parameter
 *             dependent. Internally these params are converted to data
 *             because scaling is applied.
 * @tparam T_f type of scaling vector for residual. We allow
 *             it to be @c var because scaling could be parameter
 *             dependent. Internally these params are converted to data
 *             because scaling is applied.
 * @param x initial point and final solution.
 * @param env KINSOL solution environment
 * @param f_tol Function tolerance
 * @param max_num_steps max nb. of iterations.
 */
template <typename F, typename T, typename T_u, typename T_f, typename... Args,
          require_eigen_vector_t<T>* = nullptr>
Eigen::VectorXd kinsol_solve_fp(const F& f, const T& x,
                                const double function_tolerance,
                                const double max_num_steps,
                                const std::vector<T_u>& u_scale,
                                const std::vector<T_f>& f_scale,
                                std::ostream* const msgs, const Args&... args) {
  KinsolFixedPointEnv<F, T_u, T_f, Args...> env(f, x, u_scale, f_scale, msgs,
                                                args...);
  const int N = x.size();
  constexpr int default_anderson_depth = 4;
  const int anderson_depth = std::min(N, default_anderson_depth);

  check_flag_sundials(KINSetNumMaxIters(env.mem_, max_num_steps),
                      "KINSetNumMaxIters");
  check_flag_sundials(KINSetMAA(env.mem_, anderson_depth), "KINSetMAA");
  check_flag_sundials(KINInit(env.mem_, &env.kinsol_f_system, env.nv_x_),
                      "KINInit");
  check_flag_sundials(KINSetFuncNormTol(env.mem_, function_tolerance),
                      "KINSetFuncNormTol");
  check_flag_sundials(KINSetUserData(env.mem_, static_cast<void*>(&env)),
                      "KINSetUserData");

  check_flag_kinsol(
      KINSol(env.mem_, env.nv_x_, KIN_FP, env.nv_u_scal_, env.nv_f_scal_),
      max_num_steps);

  return Eigen::Map<Eigen::VectorXd>(N_VGetArrayPointer(env.nv_x_), N);
}

/**
 * Return a fixed pointer to the specified system of algebraic
 * equations of form
 * \[
 *   x = F(x; theta)
 * \]
 * given an initial guess \(x\), and parameters \(theta\) and data. Use the
 * KINSOL solver from the SUNDIALS suite.
 *
 * The user can also specify the scaling controls, the function
 * tolerance, and the maximum number of steps.
 *
 * This function is overloaded to handle both constant and var-type parameters.
 * This overload handles var parameters, and checks the input, calls the
 * algebraic solver, and appropriately handles derivative propagation through
 * the `reverse_pass_callback`.
 *
 * @tparam F type of equation system function.
 * @tparam T type of initial guess vector.
 * @tparam T_u type of scaling vector for unknowns. We allow
 *             it to be @c var because scaling could be parameter
 *             dependent. Internally these params are converted to data
 *             because scaling is applied.
 * @tparam T_f type of scaling vector for residual. We allow
 *             it to be @c var because scaling could be parameter
 *             dependent. Internally these params are converted to data
 *             because scaling is applied.
 *
 * @param[in] f functor that evaluated the system of equations.
 * @param[in] x vector of starting values.
 * @param[in, out] msgs the print stream for warning messages.
 * @param[in] u_scale diagonal scaling matrix elements \(Du\)
 *                    such that \(Du x\) has all components roughly the same
 *                    magnitude when \(x\) is close to a solution.
 *                    (ref. KINSOL user guide chap.2 sec. "Scaling")
 * @param[in] f_scale diagonal scaling matrix elements such
 *                    that \(Df (x - f(x))\) has all components roughly the same
 *                    magnitude when \(x\) is not too close to a solution.
 *                    (ref. KINSOL user guide chap.2 sec. "Scaling")
 * @param[in] function_tolerance Function-norm stopping tolerance.
 * @param[in] max_num_steps maximum number of function evaluations.
 * @param[in] args additional parameters to the equation system functor.
 * @pre f returns finite values when passed any value of x and the given args.
 * @throw <code>std::invalid_argument</code> if x has size zero.
 * @throw <code>std::invalid_argument</code> if x has non-finite elements.
 * @throw <code>std::invalid_argument</code> if scaled_step_size is strictly
 * negative.
 * @throw <code>std::invalid_argument</code> if function_tolerance is strictly
 * negative.
 * @throw <code>std::invalid_argument</code> if max_num_steps is not positive.
 * @throw <code>boost::math::evaluation_error</code> (which is a subclass of
 * <code>std::runtime_error</code>) if solver exceeds max_num_steps.
 */
template <typename F, typename T, typename T_u, typename T_f, typename... Args,
          require_eigen_vector_t<T>* = nullptr,
          require_all_st_arithmetic<Args...>* = nullptr>
Eigen::VectorXd algebra_solver_fp_impl(const F& f, const T& x,
                                       std::ostream* const msgs,
                                       const std::vector<T_u>& u_scale,
                                       const std::vector<T_f>& f_scale,
                                       const double function_tolerance,
                                       const int max_num_steps,
                                       const Args&... args) {
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
 * Signature to maintain backward compatibility, will be removed
 * in the future.
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
 * @param[in] dat Continuous data vector for the equation system.
 * @param[in] dat_int Integer data vector for the equation system.
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
 * @pre f returns finite values when passed any value of x and the given y, dat,
 *        and dat_int.
 * @throw <code>std::invalid_argument</code> if x has size zero.
 * @throw <code>std::invalid_argument</code> if x has non-finite elements.
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
template <typename F, typename T1, typename T2, typename T_u, typename T_f,
          require_eigen_vector_t<T1>* = nullptr,
          require_eigen_vector_t<T2>* = nullptr>
Eigen::Matrix<scalar_type_t<T2>, Eigen::Dynamic, 1> algebra_solver_fp(
    const F& f, const T1& x, const T2& y, const std::vector<double>& dat,
    const std::vector<int>& dat_int, const std::vector<T_u>& u_scale,
    const std::vector<T_f>& f_scale, std::ostream* const msgs = nullptr,
    double function_tolerance = 1e-8, int max_num_steps = 200) {
  return algebra_solver_fp_impl(algebra_solver_adapter<F>(f), to_ref(x), msgs, u_scale,
                                f_scale, function_tolerance, max_num_steps, to_ref(y),
                                dat, dat_int);
}

}  // namespace math
}  // namespace stan

#endif
