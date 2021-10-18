#ifndef STAN_MATH_REV_FUNCTOR_ODE_ADJOINT_HPP
#define STAN_MATH_REV_FUNCTOR_ODE_ADJOINT_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/prim/fun/eval.hpp>
#include <stan/math/rev/functor/cvodes_integrator_adjoint.hpp>
#include <ostream>
#include <vector>

namespace stan {
namespace math {

/**
 * Solve the ODE initial value problem y' = f(t, y), y(t0) = y0 at a set of
 * times, { t1, t2, t3, ... } using the stiff backward differentiation formula
 * BDF solver or the non-stiff Adams solver from CVODES. The ODE system is
 * integrated using the adjoint sensitivity approach of CVODES.
 *
 * \p f must define an operator() with the signature as:
 *   template<typename T_t, typename T_y, typename... T_Args>
 *   Eigen::Matrix<stan::return_type_t<T_t, T_y, T_Args...>, Eigen::Dynamic, 1>
 *     operator()(const T_t& t, const Eigen::Matrix<T_y, Eigen::Dynamic, 1>& y,
 *     std::ostream* msgs, const T_Args&... args);
 *
 * t is the time, y is the state, msgs is a stream for error messages, and args
 * are optional arguments passed to the ODE solve function (which are passed
 * through to \p f without modification).
 *
 * @tparam F Type of ODE right hand side
 * @tparam T_y0 Type of initial state
 * @tparam T_t0 Type of initial time
 * @tparam T_ts Type of output times
 * @tparam T_Args Types of pass-through parameters
 *
 * @param function_name Calling function name (for printing debugging messages)
 * @param f Right hand side of the ODE
 * @param y0 Initial state
 * @param t0 Initial time
 * @param ts Times at which to solve the ODE at. All values must be sorted and
 *   not less than t0.
 * @param relative_tolerance_forward Relative tolerance for forward problem
 * passed to CVODES
 * @param absolute_tolerance_forward Absolute tolerance per ODE state for
 * forward problem passed to CVODES
 * @param relative_tolerance_backward Relative tolerance for backward problem
 * passed to CVODES
 * @param absolute_tolerance_backward Absolute tolerance per ODE state for
 * backward problem passed to CVODES
 * @param relative_tolerance_quadrature Relative tolerance for quadrature
 * problem passed to CVODES
 * @param absolute_tolerance_quadrature Absolute tolerance for quadrature
 * problem passed to CVODES
 * @param max_num_steps Upper limit on the number of integration steps to
 *   take between each output (error if exceeded)
 * @param num_steps_between_checkpoints Number of integrator steps after which a
 * checkpoint is stored for the backward pass
 * @param interpolation_polynomial type of polynomial used for interpolation
 * @param solver_forward solver used for forward pass
 * @param solver_backward solver used for backward pass
 * @param[in, out] msgs the print stream for warning messages
 * @param args Extra arguments passed unmodified through to ODE right hand side
 * @return An `std::vector` of Eigen column vectors with scalars equal to
 *  the least upper bound of `T_y0`, `T_t0`, `T_ts`, and the lambda's arguments.
 *  This represents the solution to ODE at times \p ts
 */
template <typename F, typename T_y0, typename T_t0, typename T_ts,
          typename T_abs_tol_fwd, typename T_abs_tol_bwd, typename... T_Args,
          require_all_eigen_col_vector_t<T_y0, T_abs_tol_fwd,
                                         T_abs_tol_bwd>* = nullptr,
          require_any_not_st_arithmetic<T_y0, T_t0, T_ts, T_Args...>* = nullptr>
auto ode_adjoint_impl(
    const char* function_name, F&& f, const T_y0& y0, const T_t0& t0,
    const std::vector<T_ts>& ts, double relative_tolerance_forward,
    const T_abs_tol_fwd& absolute_tolerance_forward,
    double relative_tolerance_backward,
    const T_abs_tol_bwd& absolute_tolerance_backward,
    double relative_tolerance_quadrature, double absolute_tolerance_quadrature,
    long int max_num_steps,                  // NOLINT(runtime/int)
    long int num_steps_between_checkpoints,  // NOLINT(runtime/int)
    int interpolation_polynomial, int solver_forward, int solver_backward,
    std::ostream* msgs, const T_Args&... args) {
  using integrator_vari
      = cvodes_integrator_adjoint_vari<F, plain_type_t<T_y0>, T_t0, T_ts,
                                       plain_type_t<T_Args>...>;
  auto integrator = new integrator_vari(
      function_name, std::forward<F>(f), eval(y0), t0, ts,
      relative_tolerance_forward, absolute_tolerance_forward,
      relative_tolerance_backward, absolute_tolerance_backward,
      relative_tolerance_quadrature, absolute_tolerance_quadrature,
      max_num_steps, num_steps_between_checkpoints, interpolation_polynomial,
      solver_forward, solver_backward, msgs, eval(args)...);
  return integrator->solution();
}

/**
 * Solve the ODE initial value problem y' = f(t, y), y(t0) = y0 at a set of
 * times, { t1, t2, t3, ... } using the stiff backward differentiation formula
 * BDF solver or the non-stiff Adams solver from CVODES. The ODE system is
 * integrated using the adjoint sensitivity approach of CVODES. This
 * implementation handles the case of a double return type which ensures that no
 * resources are left on the AD stack.
 *
 * \p f must define an operator() with the signature as:
 *   template<typename T_t, typename T_y, typename... T_Args>
 *   Eigen::Matrix<stan::return_type_t<T_t, T_y, T_Args...>, Eigen::Dynamic, 1>
 *     operator()(const T_t& t, const Eigen::Matrix<T_y, Eigen::Dynamic, 1>& y,
 *     std::ostream* msgs, const T_Args&... args);
 *
 * t is the time, y is the state, msgs is a stream for error messages, and args
 * are optional arguments passed to the ODE solve function (which are passed
 * through to \p f without modification).
 *
 * @tparam F Type of ODE right hand side
 * @tparam T_y0 Type of initial state
 * @tparam T_t0 Type of initial time
 * @tparam T_ts Type of output times
 * @tparam T_Args Types of pass-through parameters
 *
 * @param function_name Calling function name (for printing debugging messages)
 * @param f Right hand side of the ODE
 * @param y0 Initial state
 * @param t0 Initial time
 * @param ts Times at which to solve the ODE at. All values must be sorted and
 *   not less than t0.
 * @param relative_tolerance_forward Relative tolerance for forward problem
 * passed to CVODES
 * @param absolute_tolerance_forward Absolute tolerance per ODE state for
 * forward problem passed to CVODES
 * @param relative_tolerance_backward Relative tolerance for backward problem
 * passed to CVODES
 * @param absolute_tolerance_backward Absolute tolerance per ODE state for
 * backward problem passed to CVODES
 * @param relative_tolerance_quadrature Relative tolerance for quadrature
 * problem passed to CVODES
 * @param absolute_tolerance_quadrature Absolute tolerance for quadrature
 * problem passed to CVODES
 * @param max_num_steps Upper limit on the number of integration steps to
 *   take between each output (error if exceeded)
 * @param num_steps_between_checkpoints Number of integrator steps after which a
 * checkpoint is stored for the backward pass
 * @param interpolation_polynomial type of polynomial used for interpolation
 * @param solver_forward solver used for forward pass
 * @param solver_backward solver used for backward pass
 * @param[in, out] msgs the print stream for warning messages
 * @param args Extra arguments passed unmodified through to ODE right hand side
 * @return An `std::vector` of Eigen column vectors with scalars equal to
 *  the least upper bound of `T_y0`, `T_t0`, `T_ts`, and the lambda's arguments.
 *  This represents the solution to ODE at times \p ts
 */
template <typename F, typename T_y0, typename T_t0, typename T_ts,
          typename T_abs_tol_fwd, typename T_abs_tol_bwd, typename... T_Args,
          require_all_eigen_col_vector_t<T_y0, T_abs_tol_fwd,
                                         T_abs_tol_bwd>* = nullptr,
          require_all_st_arithmetic<T_y0, T_t0, T_ts, T_Args...>* = nullptr>
std::vector<Eigen::VectorXd> ode_adjoint_impl(
    const char* function_name, F&& f, const T_y0& y0, const T_t0& t0,
    const std::vector<T_ts>& ts, double relative_tolerance_forward,
    const T_abs_tol_fwd& absolute_tolerance_forward,
    double relative_tolerance_backward,
    const T_abs_tol_bwd& absolute_tolerance_backward,
    double relative_tolerance_quadrature, double absolute_tolerance_quadrature,
    long int max_num_steps,                  // NOLINT(runtime/int)
    long int num_steps_between_checkpoints,  // NOLINT(runtime/int)
    int interpolation_polynomial, int solver_forward, int solver_backward,
    std::ostream* msgs, const T_Args&... args) {
  std::vector<Eigen::VectorXd> ode_solution;
  {
    nested_rev_autodiff nested;

    using integrator_vari
        = cvodes_integrator_adjoint_vari<F, plain_type_t<T_y0>, T_t0, T_ts,
                                         plain_type_t<T_Args>...>;

    auto integrator = new integrator_vari(
        function_name, std::forward<F>(f), eval(y0), t0, ts,
        relative_tolerance_forward, absolute_tolerance_forward,
        relative_tolerance_backward, absolute_tolerance_backward,
        relative_tolerance_quadrature, absolute_tolerance_quadrature,
        max_num_steps, num_steps_between_checkpoints, interpolation_polynomial,
        solver_forward, solver_backward, msgs, eval(args)...);

    ode_solution = integrator->solution();
  }
  return ode_solution;
}

/**
 * Solve the ODE initial value problem y' = f(t, y), y(t0) = y0 at a set of
 * times, { t1, t2, t3, ... } using the stiff backward differentiation formula
 * BDF solver or the non-stiff Adams solver from CVODES. The ODE system is
 * integrated using the adjoint sensitivity approach of CVODES.
 *
 * \p f must define an operator() with the signature as:
 *   template<typename T_t, typename T_y, typename... T_Args>
 *   Eigen::Matrix<stan::return_type_t<T_t, T_y, T_Args...>, Eigen::Dynamic, 1>
 *     operator()(const T_t& t, const Eigen::Matrix<T_y, Eigen::Dynamic, 1>& y,
 *     std::ostream* msgs, const T_Args&... args);
 *
 * t is the time, y is the state, msgs is a stream for error messages, and args
 * are optional arguments passed to the ODE solve function (which are passed
 * through to \p f without modification).
 *
 * @tparam F Type of ODE right hand side
 * @tparam T_y0 Type of initial state
 * @tparam T_t0 Type of initial time
 * @tparam T_ts Type of output times
 * @tparam T_Args Types of pass-through parameters
 *
 * @param f Right hand side of the ODE
 * @param y0 Initial state
 * @param t0 Initial time
 * @param ts Times at which to solve the ODE at. All values must be sorted and
 *   not less than t0.
 * @param relative_tolerance_forward Relative tolerance for forward problem
 * passed to CVODES
 * @param absolute_tolerance_forward Absolute tolerance per ODE state for
 * forward problem passed to CVODES
 * @param relative_tolerance_backward Relative tolerance for backward problem
 * passed to CVODES
 * @param absolute_tolerance_backward Absolute tolerance per ODE state for
 * backward problem passed to CVODES
 * @param relative_tolerance_quadrature Relative tolerance for quadrature
 * problem passed to CVODES
 * @param absolute_tolerance_quadrature Absolute tolerance for quadrature
 * problem passed to CVODES
 * @param max_num_steps Upper limit on the number of integration steps to
 *   take between each output (error if exceeded)
 * @param num_steps_between_checkpoints Number of integrator steps after which a
 * checkpoint is stored for the backward pass
 * @param interpolation_polynomial type of polynomial used for interpolation
 * @param solver_forward solver used for forward pass
 * @param solver_backward solver used for backward pass
 * @param[in, out] msgs the print stream for warning messages
 * @param args Extra arguments passed unmodified through to ODE right hand side
 * @return An `std::vector` of Eigen column vectors with scalars equal to
 *  the least upper bound of `T_y0`, `T_t0`, `T_ts`, and the lambda's arguments.
 *  This represents the solution to ODE at times \p ts
 */
template <typename F, typename T_y0, typename T_t0, typename T_ts,
          typename T_abs_tol_fwd, typename T_abs_tol_bwd, typename... T_Args,
          require_all_eigen_col_vector_t<T_y0, T_abs_tol_fwd,
                                         T_abs_tol_bwd>* = nullptr>
auto ode_adjoint_tol_ctl(
    F&& f, const T_y0& y0, const T_t0& t0, const std::vector<T_ts>& ts,
    double relative_tolerance_forward,
    const T_abs_tol_fwd& absolute_tolerance_forward,
    double relative_tolerance_backward,
    const T_abs_tol_bwd& absolute_tolerance_backward,
    double relative_tolerance_quadrature, double absolute_tolerance_quadrature,
    long int max_num_steps,                  // NOLINT(runtime/int)
    long int num_steps_between_checkpoints,  // NOLINT(runtime/int)
    int interpolation_polynomial, int solver_forward, int solver_backward,
    std::ostream* msgs, const T_Args&... args) {
  return ode_adjoint_impl(
      "ode_adjoint_tol_ctl", std::forward<F>(f), y0, t0, ts,
      relative_tolerance_forward, absolute_tolerance_forward,
      relative_tolerance_backward, absolute_tolerance_backward,
      relative_tolerance_quadrature, absolute_tolerance_quadrature,
      max_num_steps, num_steps_between_checkpoints, interpolation_polynomial,
      solver_forward, solver_backward, msgs, args...);
}

}  // namespace math
}  // namespace stan
#endif
