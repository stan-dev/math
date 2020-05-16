#ifndef STAN_MATH_PRIM_FUNCTOR_ODE_RK45_HPP
#define STAN_MATH_PRIM_FUNCTOR_ODE_RK45_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/functor/coupled_ode_system.hpp>
#include <stan/math/prim/functor/ode_store_sensitivities.hpp>
#include <boost/numeric/odeint/external/eigen/eigen_algebra.hpp>
#include <boost/numeric/odeint.hpp>
#include <ostream>
#include <vector>

namespace stan {
namespace math {

/**
 * Solve the ODE initial value problem y' = f(t, y), y(t0) = y0 at a set of
 * times, { t1, t2, t3, ... } using the non-stiff Runge-Kutta 45 solver in Boost
 * with a relative tolerance of 1e-10, an absolute tolerance of 1e-10, and
 * taking a maximum of 1e8 steps.
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
template <typename F, typename T_initial, typename T_t0, typename T_ts,
          typename... Args>
std::vector<
    Eigen::Matrix<stan::return_type_t<T_initial, Args...>, Eigen::Dynamic, 1>>
ode_rk45(const F& f, const Eigen::Matrix<T_initial, Eigen::Dynamic, 1>& y0,
         double t0, const std::vector<double>& ts, std::ostream* msgs,
         const Args&... args) {
  double relative_tolerance = 1e-10;
  double absolute_tolerance = 1e-10;
  long int max_num_steps = 1e8;

  return ode_rk45_tol(f, y0, t0, ts, relative_tolerance, absolute_tolerance,
                      max_num_steps, msgs, args...);
}

/**
 * Solve the ODE initial value problem y' = f(t, y), y(t0) = y0 at a set of
 * times, { t1, t2, t3, ... } using the non-stiff Runge-Kutta 45 solver in Boost
 * with a relative tolerance of 1e-10, an absolute tolerance of 1e-10, and
 * taking a maximum of 1e8 steps.
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
template <typename F, typename T_initial, typename T_t0, typename T_ts,
          typename... Args>
std::vector<Eigen::Matrix<stan::return_type_t<T_initial, T_t0, T_ts, Args...>,
                          Eigen::Dynamic, 1>>
ode_rk45_tol(const F& f,
             const Eigen::Matrix<T_initial, Eigen::Dynamic, 1>& y0_arg, T_t0 t0,
             const std::vector<T_ts>& ts, double relative_tolerance,
             double absolute_tolerance, long int max_num_steps,
             std::ostream* msgs, const Args&... args) {
  using boost::numeric::odeint::integrate_times;
  using boost::numeric::odeint::make_dense_output;
  using boost::numeric::odeint::max_step_checker;
  using boost::numeric::odeint::runge_kutta_dopri5;
  using boost::numeric::odeint::vector_space_algebra;

  using T_initial_or_t0 = return_type_t<T_initial, T_t0>;

  Eigen::Matrix<T_initial_or_t0, Eigen::Dynamic, 1> y0 = y0_arg.unaryExpr(
      [](const T_initial& val) { return T_initial_or_t0(val); });
  const std::vector<double> ts_dbl = value_of(ts);

  check_finite("integrate_ode_rk45", "initial state", y0);
  check_finite("integrate_ode_rk45", "initial time", t0);
  check_finite("integrate_ode_rk45", "times", ts_dbl);

  // Code from https://stackoverflow.com/a/17340003
  std::vector<int> unused_temp{
      0, (check_finite("integrate_ode_rk45", "ode parameters and data", args),
          0)...};

  check_nonzero_size("integrate_ode_rk45", "initial state", y0);
  check_nonzero_size("integrate_ode_rk45", "times", ts_dbl);
  check_ordered("integrate_ode_rk45", "times", ts_dbl);
  check_less("integrate_ode_rk45", "initial time", t0, ts_dbl[0]);

  if (relative_tolerance <= 0) {
    invalid_argument("integrate_ode_rk45", "relative_tolerance,",
                     relative_tolerance, "", ", must be greater than 0");
  }
  if (absolute_tolerance <= 0) {
    invalid_argument("integrate_ode_rk45", "absolute_tolerance,",
                     absolute_tolerance, "", ", must be greater than 0");
  }
  if (max_num_steps <= 0) {
    invalid_argument("integrate_ode_rk45", "max_num_steps,", max_num_steps, "",
                     ", must be greater than 0");
  }

  using return_t = return_type_t<T_initial, T_t0, T_ts, Args...>;
  // creates basic or coupled system by template specializations
  coupled_ode_system<F, T_initial_or_t0, Args...> coupled_system(f, y0, msgs,
                                                                 args...);

  // first time in the vector must be time of initial state
  std::vector<double> ts_vec(ts.size() + 1);
  ts_vec[0] = value_of(t0);
  std::copy(ts_dbl.begin(), ts_dbl.end(), ts_vec.begin() + 1);

  std::vector<Eigen::Matrix<return_t, Eigen::Dynamic, 1>> y;
  bool observer_initial_recorded = false;
  size_t time_index = 0;

  // avoid recording of the initial state which is included by the
  // conventions of odeint in the output
  auto filtered_observer
      = [&](const Eigen::VectorXd& coupled_state, double t) -> void {
    if (!observer_initial_recorded) {
      observer_initial_recorded = true;
      return;
    }
    y.emplace_back(ode_store_sensitivities(f, coupled_state, y0, t0,
                                           ts[time_index], msgs, args...));
    time_index++;
  };

  // the coupled system creates the coupled initial state
  Eigen::VectorXd initial_coupled_state = coupled_system.initial_state();

  const double step_size = 0.1;
  integrate_times(
      make_dense_output(
          absolute_tolerance, relative_tolerance,
          runge_kutta_dopri5<Eigen::VectorXd, double, Eigen::VectorXd, double,
                             vector_space_algebra>()),
      std::ref(coupled_system), initial_coupled_state, std::begin(ts_vec),
      std::end(ts_vec), step_size, filtered_observer,
      max_step_checker(max_num_steps));

  return y;
}

}  // namespace math
}  // namespace stan
#endif
