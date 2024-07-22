#ifndef STAN_MATH_PRIM_FUNCTOR_ODE_CKRK_HPP
#define STAN_MATH_PRIM_FUNCTOR_ODE_CKRK_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/functor/apply.hpp>
#include <stan/math/prim/functor/coupled_ode_system.hpp>
#include <stan/math/prim/functor/ode_store_sensitivities.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <boost/numeric/odeint.hpp>
#include <ostream>
#include <tuple>
#include <vector>

namespace stan {
namespace math {

/**
 * Solve the ODE initial value problem y' = f(t, y), y(t0) = y0 at a set of
 * times, { t1, t2, t3, ... } using Boost's Cash-Karp54 solver.
 *
 * If the system of equations is stiff, <code>ode_bdf</code> will likely be
 * faster.
 *
 * \p f must define an operator() with the signature as:
 *   template<typename T_t, typename T_y, typename... T_Args>
 *   Eigen::Matrix<stan::return_type_t<T_t, T_y, T_Args...>, Eigen::Dynamic, 1>
 *     operator()(const T_t& t, const Eigen::Matrix<T_y, Eigen::Dynamic, 1>& y,
 *     std::ostream* msgs, const T_Args&... args);
 *
 * t is the time, y is the vector-valued state, msgs is a stream for error
 * messages, and args are optional arguments passed to the ODE solve function
 * (which are passed through to \p f without modification).
 *
 * @tparam F Type of ODE right hand side
 * @tparam T_y0 Type of initial condition
 * @tparam T_t0 Type of initial time
 * @tparam T_ts Type of output times
 * @tparam T_Args Types of pass-through parameters
 *
 * @param function_name Calling function name (for printing debugging messages)
 * @param f Right hand side of the ODE
 * @param y0_arg Initial state
 * @param t0 Initial time
 * @param ts Times at which to solve the ODE at. All values must be sorted and
 *   greater than t0.
 * @param relative_tolerance Relative tolerance passed to Boost
 * @param absolute_tolerance Absolute tolerance passed to Boost
 * @param max_num_steps Upper limit on the number of integration steps to
 *   take between each output (error if exceeded)
 * @param[in, out] msgs the print stream for warning messages
 * @param args Extra arguments passed unmodified through to ODE right hand side
 * @return Solution to ODE at times \p ts
 */
template <typename F, typename T_y0, typename T_t0, typename T_ts,
          typename... Args, require_eigen_vector_t<T_y0>* = nullptr>
std::vector<Eigen::Matrix<stan::return_type_t<T_y0, T_t0, T_ts, Args...>,
                          Eigen::Dynamic, 1>>
ode_ckrk_tol_impl(const char* function_name, const F& f, const T_y0& y0_arg,
                  T_t0 t0, const std::vector<T_ts>& ts,
                  double relative_tolerance, double absolute_tolerance,
                  long int max_num_steps,  // NOLINT(runtime/int)
                  std::ostream* msgs, const Args&... args) {
  using boost::numeric::odeint::integrate_times;
  using boost::numeric::odeint::make_dense_output;
  using boost::numeric::odeint::max_step_checker;
  using boost::numeric::odeint::no_progress_error;
  using boost::numeric::odeint::runge_kutta_cash_karp54;
  using boost::numeric::odeint::vector_space_algebra;

  using T_y0_t0 = return_type_t<T_y0, T_t0>;

  Eigen::Matrix<T_y0_t0, Eigen::Dynamic, 1> y0
      = y0_arg.template cast<T_y0_t0>();

  check_finite(function_name, "initial state", y0);
  check_finite(function_name, "initial time", t0);
  check_finite(function_name, "times", ts);

  std::tuple<ref_type_t<Args>...> args_ref_tuple(args...);

  math::apply(
      [&](const auto&... args_ref) {
        // Code from https://stackoverflow.com/a/17340003
        std::vector<int> unused_temp{
            0,
            (check_finite(function_name, "ode parameters and data", args_ref),
             0)...};
      },
      args_ref_tuple);

  check_nonzero_size(function_name, "initial state", y0);
  check_nonzero_size(function_name, "times", ts);
  check_sorted(function_name, "times", ts);
  check_less(function_name, "initial time", t0, ts[0]);

  check_positive_finite(function_name, "relative_tolerance",
                        relative_tolerance);
  check_positive_finite(function_name, "absolute_tolerance",
                        absolute_tolerance);
  check_positive(function_name, "max_num_steps", max_num_steps);

  using return_t = return_type_t<T_y0, T_t0, T_ts, Args...>;
  // creates basic or coupled system by template specializations
  auto&& coupled_system = math::apply(
      [&](const auto&... args_ref) {
        return coupled_ode_system<F, T_y0_t0, ref_type_t<Args>...>(f, y0, msgs,
                                                                   args_ref...);
      },
      args_ref_tuple);

  // first time in the vector must be time of initial state
  std::vector<double> ts_vec(ts.size() + 1);
  ts_vec[0] = value_of(t0);
  for (size_t i = 0; i < ts.size(); ++i)
    ts_vec[i + 1] = value_of(ts[i]);

  std::vector<Eigen::Matrix<return_t, Eigen::Dynamic, 1>> y;
  y.reserve(ts.size());
  bool observer_initial_recorded = false;
  size_t time_index = 0;

  // avoid recording of the initial state which is included by the
  // conventions of odeint in the output
  auto filtered_observer
      = [&](const std::vector<double>& coupled_state, double t) -> void {
    if (!observer_initial_recorded) {
      observer_initial_recorded = true;
      return;
    }
    math::apply(
        [&](const auto&... args_ref) {
          y.emplace_back(ode_store_sensitivities(
              f, coupled_state, y0, t0, ts[time_index], msgs, args_ref...));
        },
        args_ref_tuple);
    time_index++;
  };

  // the coupled system creates the coupled initial state
  std::vector<double> initial_coupled_state = coupled_system.initial_state();

  const double step_size = 0.1;
  try {
    integrate_times(
        make_controlled(absolute_tolerance, relative_tolerance,
                        runge_kutta_cash_karp54<std::vector<double>, double,
                                                std::vector<double>, double>()),
        std::ref(coupled_system), initial_coupled_state, std::begin(ts_vec),
        std::end(ts_vec), step_size, filtered_observer,
        max_step_checker(max_num_steps));
  } catch (const no_progress_error& e) {
    throw_domain_error(function_name, "", ts_vec[time_index + 1],
                       "Failed to integrate to next output time (",
                       ") in less than max_num_steps steps");
  }

  return y;
}

/**
 * Solve the ODE initial value problem y' = f(t, y), y(t0) = y0 at a set of
 * times, { t1, t2, t3, ... } using Boost's Cash-Karp solver.
 *
 * If the system of equations is stiff, <code>ode_bdf</code> will likely be
 * faster.
 *
 * \p f must define an operator() with the signature as:
 *   template<typename T_t, typename T_y, typename... T_Args>
 *   Eigen::Matrix<stan::return_type_t<T_t, T_y, T_Args...>, Eigen::Dynamic, 1>
 *     operator()(const T_t& t, const Eigen::Matrix<T_y, Eigen::Dynamic, 1>& y,
 *     std::ostream* msgs, const T_Args&... args);
 *
 * t is the time, y is the vector-valued state, msgs is a stream for error
 * messages, and args are optional arguments passed to the ODE solve function
 * (which are passed through to \p f without modification).
 *
 * @tparam F Type of ODE right hand side
 * @tparam T_0 Type of initial time
 * @tparam T_ts Type of output times
 * @tparam T_Args Types of pass-through parameters
 *
 * @param f Right hand side of the ODE
 * @param y0_arg Initial state
 * @param t0 Initial time
 * @param ts Times at which to solve the ODE at. All values must be sorted and
 *   greater than t0.
 * @param relative_tolerance Relative tolerance passed to Boost
 * @param absolute_tolerance Absolute tolerance passed to Boost
 * @param max_num_steps Upper limit on the number of integration steps to
 *   take between each output (error if exceeded)
 * @param[in, out] msgs the print stream for warning messages
 * @param args Extra arguments passed unmodified through to ODE right hand side
 * @return Solution to ODE at times \p ts
 */
template <typename F, typename T_y0, typename T_t0, typename T_ts,
          typename... Args, require_eigen_vector_t<T_y0>* = nullptr>
std::vector<Eigen::Matrix<stan::return_type_t<T_y0, T_t0, T_ts, Args...>,
                          Eigen::Dynamic, 1>>
ode_ckrk_tol(const F& f, const T_y0& y0_arg, T_t0 t0,
             const std::vector<T_ts>& ts, double relative_tolerance,
             double absolute_tolerance,
             long int max_num_steps,  // NOLINT(runtime/int)
             std::ostream* msgs, const Args&... args) {
  return ode_ckrk_tol_impl("ode_ckrk_tol", f, y0_arg, t0, ts,
                           relative_tolerance, absolute_tolerance,
                           max_num_steps, msgs, args...);
}

/**
 * Solve the ODE initial value problem y' = f(t, y), y(t0) = y0 at a set of
 * times, { t1, t2, t3, ... } using Boost's Cash-Karp Runge-Kutta solver
 * with defaults for relative_tolerance, absolute_tolerance, and max_num_steps.
 *
 * If the system of equations is stiff, <code>ode_bdf</code> will likely be
 * faster.
 *
 * \p f must define an operator() with the signature as:
 *   template<typename T_t, typename T_y, typename... T_Args>
 *   Eigen::Matrix<stan::return_type_t<T_t, T_y, T_Args...>, Eigen::Dynamic, 1>
 *     operator()(const T_t& t, const Eigen::Matrix<T_y, Eigen::Dynamic, 1>& y,
 *     std::ostream* msgs, const T_Args&... args);
 *
 * t is the time, y is the vector-valued state, msgs is a stream for error
 * messages, and args are optional arguments passed to the ODE solve function
 * (which are passed through to \p f without modification).
 *
 * @tparam F Type of ODE right hand side
 * @tparam T_y0 Type of initial time
 * @tparam T_ts Type of output times
 * @tparam T_Args Types of pass-through parameters
 *
 * @param f Right hand side of the ODE
 * @param y0 Initial state
 * @param t0 Initial time
 * @param ts Times at which to solve the ODE at. All values must be sorted and
 *   greather than t0.
 * @param[in, out] msgs the print stream for warning messages
 * @param args Extra arguments passed unmodified through to ODE right hand side
 * @return Solution to ODE at times \p ts
 */
template <typename F, typename T_y0, typename T_t0, typename T_ts,
          typename... Args, require_eigen_vector_t<T_y0>* = nullptr>
std::vector<Eigen::Matrix<stan::return_type_t<T_y0, T_t0, T_ts, Args...>,
                          Eigen::Dynamic, 1>>
ode_ckrk(const F& f, const T_y0& y0, T_t0 t0, const std::vector<T_ts>& ts,
         std::ostream* msgs, const Args&... args) {
  double relative_tolerance = 1e-6;
  double absolute_tolerance = 1e-6;
  long int max_num_steps = 1e6;  // NOLINT(runtime/int)

  return ode_ckrk_tol_impl("ode_ckrk", f, y0, t0, ts, relative_tolerance,
                           absolute_tolerance, max_num_steps, msgs, args...);
}

}  // namespace math
}  // namespace stan
#endif
