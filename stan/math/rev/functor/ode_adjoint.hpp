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
 * BDF solver from CVODES.
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
 * @tparam T_0 Type of initial time
 * @tparam T_ts Type of output times
 * @tparam T_Args Types of pass-through parameters
 *
 * @param function_name Calling function name (for printing debugging messages)
 * @param f Right hand side of the ODE
 * @param y0 Initial state
 * @param t0 Initial time
 * @param ts Times at which to solve the ODE at. All values must be sorted and
 *   not less than t0.
 * @param rel_tol_f Relative tolerance for forward problem passed to CVODES
 * @param abs_tol_f Absolute tolerance for forward problem passed to CVODES
 * @param rel_tol_b Relative tolerance for backward problem passed to CVODES
 * @param abs_tol_b Absolute tolerance for backward problem passed to CVODES
 * @param rel_tol_q Relative tolerance for quadrature problem passed to CVODES
 * @param abs_tol_q Absolute tolerance for quadrature problem passed to CVODES
 * @param max_num_steps Upper limit on the number of integration steps to
 *   take between each output (error if exceeded)
 * @param num_checkpoints Number of integrator steps after which a checkpoint is stored for the backward pass
 * @param interpolation_polynomial type of polynomial used for interpolation
 * @param solver_f solver used for forward pass
 * @param solver_b solver used for backward pass
 * @param[in, out] msgs the print stream for warning messages
 * @param args Extra arguments passed unmodified through to ODE right hand side
 * @return Solution to ODE at times \p ts
 */
template <typename F, typename T_y0, typename T_t0, typename T_ts,
          typename... T_Args, require_eigen_col_vector_t<T_y0>* = nullptr>
std::vector<Eigen::Matrix<stan::return_type_t<T_y0, T_t0, T_ts, T_Args...>,
                          Eigen::Dynamic, 1>>
ode_adjoint_impl(const char* function_name, const F& f,
                 const T_y0& y0,
                 const T_t0& t0, const std::vector<T_ts>& ts,
                 double rel_tol_f, Eigen::VectorXd abs_tol_f,
                 double rel_tol_b, double abs_tol_b,
                 double rel_tol_q, double abs_tol_q,
                 long int max_num_steps, long int num_checkpoints,
                 int interpolation_polynomial,
                 int solver_f, int solver_b,
                 std::ostream* msgs,
                 const T_Args&... args) {
  /*
  const auto& args_ref_tuple = std::make_tuple(to_ref(args)...);
  return apply(
      [&](const auto&... args_refs) {
        auto integrator
            = new stan::math::cvodes_integrator_adjoint_vari<CV_BDF, F, T_y0,
  T_t0, T_ts, T_Args...>( function_name, f, y0, t0, ts, relative_tolerance,
                          absolute_tolerance, absolute_tolerance_B,
  absolute_tolerance_QB, steps_checkpoint, max_num_steps, msgs, args_refs...);
        return (*integrator)();
      }, args_ref_tuple);
      */
  /*
  const auto& args_eval_tuple = std::make_tuple(eval(args)...);
  return apply(
      [&](const auto&... args_eval) {
        auto integrator
            = new stan::math::cvodes_integrator_adjoint_vari<CV_BDF, F, T_y0,
  T_t0, T_ts, T_Args...>( function_name, f, y0, t0, ts, relative_tolerance,
                          absolute_tolerance, absolute_tolerance_B,
  absolute_tolerance_QB, steps_checkpoint, max_num_steps, msgs, args_eval...);
        return (*integrator)();
      }, args_eval_tuple);
      */
  /*
  auto integrator
      = new stan::math::cvodes_integrator_adjoint_vari<CV_BDF, F, T_y0, T_t0,
                                                       T_ts, T_Args...>(
          function_name, f, y0, t0, ts, relative_tolerance, absolute_tolerance,
          absolute_tolerance_B, absolute_tolerance_QB, steps_checkpoint,
          max_num_steps, msgs, eval(args)...);
  return (*integrator)();
  */
  auto integrator = new stan::math::cvodes_integrator_adjoint_vari<
      F, plain_type_t<T_y0>, T_t0, T_ts, plain_type_t<T_Args>...>(
      function_name, f, eval(y0), t0, ts,
      rel_tol_f, abs_tol_f,
      rel_tol_b, abs_tol_b,
      rel_tol_q, abs_tol_q,
      max_num_steps, num_checkpoints, interpolation_polynomial, solver_f, solver_b,
      msgs,
      eval(args)...);
  return (*integrator)();
}

/**
 * Solve the ODE initial value problem y' = f(t, y), y(t0) = y0 at a set of
 * times, { t1, t2, t3, ... } using the stiff backward differentiation formula
 * BDF solver from CVODES.
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
 * @tparam T_0 Type of initial time
 * @tparam T_ts Type of output times
 * @tparam T_Args Types of pass-through parameters
 *
 * @param f Right hand side of the ODE
 * @param y0 Initial state
 * @param t0 Initial time
 * @param ts Times at which to solve the ODE at. All values must be sorted and
 *   not less than t0.
 * @param relative_tolerance Relative tolerance passed to CVODES
 * @param absolute_tolerance Absolute tolerance passed to CVODES
 * @param max_num_steps Upper limit on the number of integration steps to
 *   take between each output (error if exceeded)
 * @param[in, out] msgs the print stream for warning messages
 * @param args Extra arguments passed unmodified through to ODE right hand side
 * @return Solution to ODE at times \p ts
 */
template <typename F, typename T_y0, typename T_t0, typename T_ts,
          typename... T_Args, require_eigen_col_vector_t<T_y0>* = nullptr>
std::vector<Eigen::Matrix<stan::return_type_t<T_y0, T_t0, T_ts, T_Args...>,
                          Eigen::Dynamic, 1>>
ode_adjoint(const F& f, const T_y0& y0, const T_t0& t0,
            const std::vector<T_ts>& ts,
            std::ostream* msgs, // @rok: should be moved to end like
                                // for bdf??
            double rel_tol_f, Eigen::VectorXd abs_tol_f,
            double rel_tol_b, double abs_tol_b,
            double rel_tol_q, double abs_tol_q,
            long int max_num_steps, long int num_checkpoints,
            int interpolation_polynomial,
            int solver_f, int solver_b,
            //std::ostream* msgs,
            const T_Args&... args
            ) {
  return ode_adjoint_impl(
      "ode_adjoint",
      f, y0, t0, ts,
      rel_tol_f, abs_tol_f,
      rel_tol_b, abs_tol_b,
      rel_tol_q, abs_tol_q,
      max_num_steps, num_checkpoints, interpolation_polynomial, solver_f, solver_b,
      msgs,
      args...);
}

// TODO(wds15): ask @rok/all to name this ode_adjoint_tol ?

}  // namespace math
}  // namespace stan
#endif
