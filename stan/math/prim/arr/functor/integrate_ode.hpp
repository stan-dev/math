#ifndef STAN_MATH_PRIM_ARR_FUNCTOR_INTEGRATE_ODE_HPP
#define STAN_MATH_PRIM_ARR_FUNCTOR_INTEGRATE_ODE_HPP

#include <stan/math/prim/scal/fun/value_of.hpp>
#include <stan/math/prim/scal/err/check_less.hpp>
#include <stan/math/prim/scal/err/check_finite.hpp>
#include <stan/math/prim/scal/err/check_nonzero_size.hpp>
#include <stan/math/prim/mat/err/check_ordered.hpp>
#include <stan/math/prim/scal/meta/return_type.hpp>
#include <stan/math/prim/arr/functor/coupled_ode_system.hpp>
#include <stan/math/prim/arr/functor/coupled_ode_observer.hpp>
#include <boost/numeric/odeint.hpp>
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
     * @param[in] f functor for the base ordinary differential equation.
     * @param[in] y0 initial state.
     * @param[in] t0 initial time.
     * @param[in] ts times of the desired solutions, in strictly
     * increasing order, all greater than the initial time.
     * @param[in] theta parameter vector for the ODE.
     * @param[in] x continuous data vector for the ODE.
     * @param[in] x_int integer data vector for the ODE.
     * @param[in, out] msgs the print stream for warning messages.
     * @return a vector of states, each state being a vector of the
     * same size as the state variable, corresponding to a time in ts.
     */
    template <typename F, typename T1, typename T2>
    std::vector<std::vector<typename stan::return_type<T1, T2>::type> >
    integrate_ode(const F& f,
                  const std::vector<T1> y0,
                  const double t0,
                  const std::vector<double>& ts,
                  const std::vector<T2>& theta,
                  const std::vector<double>& x,
                  const std::vector<int>& x_int,
                  std::ostream* msgs) {
      using boost::numeric::odeint::integrate_times;
      using boost::numeric::odeint::runge_kutta_dopri5;
      using boost::numeric::odeint::controlled_runge_kutta;
      using boost::numeric::odeint::make_controlled;

      stan::math::check_finite("integrate_ode", "initial state", y0);
      stan::math::check_finite("integrate_ode", "initial time", t0);
      stan::math::check_finite("integrate_ode", "times", ts);
      stan::math::check_finite("integrate_ode", "parameter vector", theta);
      stan::math::check_finite("integrate_ode", "continuous data", x);

      stan::math::check_nonzero_size("integrate_ode", "times", ts);
      stan::math::check_nonzero_size("integrate_ode", "initial state", y0);
      stan::math::check_ordered("integrate_ode", "times", ts);
      stan::math::check_less("integrate_ode", "initial time", t0, ts[0]);

      const double absolute_tolerance = 1e-6;
      const double relative_tolerance = 1e-6;
      const double step_size = 0.1;

      // creates basic or coupled system by template specializations
      coupled_ode_system<F, T1, T2>
        coupled_system(f, y0, theta, x, x_int, msgs);

      std::vector<std::vector<double> > y_coupled(ts.size());

      // the coupled system creates the coupled initial state
      std::vector<double> coupled_state
        = coupled_system.initial_state();

      typedef runge_kutta_dopri5<std::vector<double>, double, std::vector<double>, double> stepper_type;
      typedef controlled_runge_kutta< stepper_type > controlled_stepper_type;

      controlled_stepper_type rk45 = make_controlled<stepper_type>(absolute_tolerance, relative_tolerance);

      integrate_adaptive(rk45, coupled_system, coupled_state, t0, ts[0], step_size);
      y_coupled[0] = coupled_state;
      
      for (size_t t = 0; t < ts.size() - 1; t++) {
	integrate_adaptive(rk45, coupled_system, coupled_state, ts[t], ts[t+1], step_size);
	y_coupled[t+1] = coupled_state;
      }

      // the coupled system also encapsulates the decoupling operation
      return coupled_system.decouple_states(y_coupled);
    }

  }

}

#endif
