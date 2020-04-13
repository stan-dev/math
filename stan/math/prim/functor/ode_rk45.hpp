#ifndef STAN_MATH_PRIM_FUNCTOR_ODE_RK45_HPP
#define STAN_MATH_PRIM_FUNCTOR_ODE_RK45_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/functor/coupled_ode_system.hpp>
#include <stan/math/prim/fun/ode_store_sensitivities.hpp>
#include <boost/numeric/odeint.hpp>
#include <ostream>
#include <vector>

namespace stan {
namespace math {

template <typename F, typename T_initial, typename... Args>
std::vector<std::vector<typename stan::return_type<T_initial, Args...>::type>>
ode_rk45(const F& f, const std::vector<T_initial>& y0,
	 double t0,
	 const std::vector<double>& ts,
	 std::ostream* msgs,
	 const Args&... args) {
  double relative_tolerance = 1e-6;
  double absolute_tolerance = 1e-6;
  long int max_num_steps = 1e6;

  return ode_rk45_tol(f, y0, t0, ts,
		      relative_tolerance, absolute_tolerance,
		      max_num_steps,
		      msgs,
		      args...);
}

template <typename F, typename T_initial, typename... Args>
std::vector<std::vector<typename stan::return_type<T_initial, Args...>::type>>
ode_rk45_tol(const F& f, const std::vector<T_initial>& y0,
	     double t0,
	     const std::vector<double>& ts,
	     double relative_tolerance,
	     double absolute_tolerance,
	     long int max_num_steps,
	     std::ostream* msgs,
	     const Args&... args) {
  using boost::numeric::odeint::integrate_times;
  using boost::numeric::odeint::make_dense_output;
  using boost::numeric::odeint::max_step_checker;
  using boost::numeric::odeint::runge_kutta_dopri5;

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

  using return_t = return_type_t<T_initial, Args...>;
  // creates basic or coupled system by template specializations
  coupled_ode_system<F, T_initial, Args...> coupled_system(f, y0, msgs, args...);

  // first time in the vector must be time of initial state
  std::vector<double> ts_vec(ts.size() + 1);
  ts_vec[0] = t0;
  std::copy(ts_dbl.begin(), ts_dbl.end(), ts_vec.begin() + 1);

  std::vector<std::vector<return_t>> y;
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
    y.emplace_back(ode_store_sensitivities(coupled_state, y0, args...));
    time_index++;
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
