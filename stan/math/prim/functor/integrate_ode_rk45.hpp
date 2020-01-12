#ifndef STAN_MATH_PRIM_FUNCTOR_INTEGRATE_ODE_RK45_HPP
#define STAN_MATH_PRIM_FUNCTOR_INTEGRATE_ODE_RK45_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/functor/coupled_ode_system.hpp>
#include <stan/math/prim/functor/coupled_ode_observer.hpp>
#include <stan/math/prim/fun/value_of.hpp>
#include <boost/version.hpp>
#if BOOST_VERSION == 106400
#include <boost/serialization/array_wrapper.hpp>
#endif
#include <boost/numeric/odeint.hpp>
#include <algorithm>
#include <ostream>
#include <functional>
#include <iterator>
#include <vector>

namespace stan {
namespace math {

/**
 * Return the solutions for the specified system of ordinary
 * differential equations given the specified initial state,
 * initial times, times of desired solution, and parameters and
 * data, writing error and warning messages to the specified
 * stream.
 *
 * <b>Warning:</b> If the system of equations is stiff, roughly
 * defined by having varying time scales across dimensions, then
 * this solver is likely to be slow.
 *
 * This function is templated to allow the initial times to be
 * either data or autodiff variables and the parameters to be data
 * or autodiff variables.  The autodiff-based implementation for
 * reverse-mode are defined in namespace <code>stan::math</code>
 * and may be invoked via argument-dependent lookup by including
 * their headers.
 *
 * This function uses the <a
 * href="http://en.wikipedia.org/wiki/Dormand–Prince_method">Dormand-Prince
 * method</a> as implemented in Boost's <code>
 * boost::numeric::odeint::runge_kutta_dopri5</code> integrator.
 *
 * @tparam F type of ODE system function.
 * @tparam T1 type of scalars for initial values.
 * @tparam T2 type of scalars for parameters.
 * @tparam T_t0 type of scalar of initial time point.
 * @tparam T_ts type of time-points where ODE solution is returned.
 * @param[in] f functor for the base ordinary differential equation.
 * @param[in] y0 initial state.
 * @param[in] t0 initial time.
 * @param[in] ts times of the desired solutions, in strictly
 * increasing order, all greater than the initial time.
 * @param[in] theta parameter vector for the ODE.
 * @param[in] x continuous data vector for the ODE.
 * @param[in] x_int integer data vector for the ODE.
 * @param[out] msgs the print stream for warning messages.
 * @param[in] relative_tolerance relative tolerance parameter
 *   for Boost's ode solver. Defaults to 1e-6.
 * @param[in] absolute_tolerance absolute tolerance parameter
 *   for Boost's ode solver. Defaults to 1e-6.
 * @param[in] max_num_steps maximum number of steps to take within
 *   the Boost ode solver.
 * @return a vector of states, each state being a vector of the
 * same size as the state variable, corresponding to a time in ts.
 */
template <typename F, typename T1, typename T2, typename T_t0, typename T_ts>
std::vector<std::vector<return_type_t<T1, T2, T_t0, T_ts>>> integrate_ode_rk45(
    const F& f, const std::vector<T1>& y0, const T_t0& t0,
    const std::vector<T_ts>& ts, const std::vector<T2>& theta,
    const std::vector<double>& x, const std::vector<int>& x_int,
    std::ostream* msgs = nullptr, double relative_tolerance = 1e-6,
    double absolute_tolerance = 1e-6, int max_num_steps = 1E6) {
  using boost::numeric::odeint::integrate_times;
  using boost::numeric::odeint::make_dense_output;
  using boost::numeric::odeint::max_step_checker;
  using boost::numeric::odeint::runge_kutta_dopri5;

  const double t0_dbl = value_of(t0);
  const std::vector<double> ts_dbl = value_of(ts);

  check_finite("integrate_ode_rk45", "initial state", y0);
  check_finite("integrate_ode_rk45", "initial time", t0_dbl);
  check_finite("integrate_ode_rk45", "times", ts_dbl);
  check_finite("integrate_ode_rk45", "parameter vector", theta);
  check_finite("integrate_ode_rk45", "continuous data", x);

  check_nonzero_size("integrate_ode_rk45", "initial state", y0);
  check_nonzero_size("integrate_ode_rk45", "times", ts_dbl);
  check_ordered("integrate_ode_rk45", "times", ts_dbl);
  check_less("integrate_ode_rk45", "initial time", t0_dbl, ts_dbl[0]);

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

  // creates basic or coupled system by template specializations
  coupled_ode_system<F, T1, T2> coupled_system(f, y0, theta, x, x_int, msgs);

  // first time in the vector must be time of initial state
  std::vector<double> ts_vec(ts.size() + 1);
  ts_vec[0] = t0_dbl;
  std::copy(ts_dbl.begin(), ts_dbl.end(), ts_vec.begin() + 1);

  std::vector<std::vector<return_type_t<T1, T2, T_t0, T_ts>>> y;
  coupled_ode_observer<F, T1, T2, T_t0, T_ts> observer(f, y0, theta, t0, ts, x,
                                                       x_int, msgs, y);
  bool observer_initial_recorded = false;

  // avoid recording of the initial state which is included by the
  // conventions of odeint in the output
  auto filtered_observer
      = [&](const std::vector<double>& coupled_state, double t) -> void {
    if (!observer_initial_recorded) {
      observer_initial_recorded = true;
      return;
    }
    observer(coupled_state, t);
  };

  // the coupled system creates the coupled initial state
  std::vector<double> initial_coupled_state = coupled_system.initial_state();

  const double step_size = 0.1;
  integrate_times(
      make_dense_output(absolute_tolerance, relative_tolerance,
                        runge_kutta_dopri5<std::vector<double>, double,
                                           std::vector<double>, double>()),
      std::ref(coupled_system), initial_coupled_state, std::begin(ts_vec),
      std::end(ts_vec), step_size, filtered_observer,
      max_step_checker(max_num_steps));

  return y;
}

}  // namespace math

}  // namespace stan

#endif
