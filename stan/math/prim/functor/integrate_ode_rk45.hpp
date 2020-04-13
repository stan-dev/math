#ifndef STAN_MATH_PRIM_FUNCTOR_INTEGRATE_ODE_RK45_HPP
#define STAN_MATH_PRIM_FUNCTOR_INTEGRATE_ODE_RK45_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/functor/ode_rk45.hpp>
#include <ostream>
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

template <typename F, typename T1, typename T_param, typename T_t0, typename T_ts>
std::vector<std::vector<return_type_t<T1, T_t0, T_ts, T_param>>> integrate_ode_rk45(
    const F& f, const std::vector<T1>& y0, const T_t0& t0,
    const std::vector<T_ts>& ts,
    const std::vector<T_param>& theta,
    const std::vector<double>& x, const std::vector<int>& x_int,
    std::ostream* msgs = nullptr,
    double relative_tolerance = 1e-6,
    double absolute_tolerance = 1e-6,
    int max_num_steps = 1e6) {
  return ode_rk45_tol(f, y0, t0, ts,
		      relative_tolerance, absolute_tolerance,
		      max_num_steps,
		      msgs,
		      theta, x, x_int);
}

}  // namespace math

}  // namespace stan

#endif
