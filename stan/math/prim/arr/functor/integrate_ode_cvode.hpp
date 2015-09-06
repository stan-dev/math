#ifndef STAN_MATH_PRIM_ARR_FUNCTOR_INTEGRATE_ODE_CVODE_HPP
#define STAN_MATH_PRIM_ARR_FUNCTOR_INTEGRATE_ODE_CVODE_HPP

#include <stan/math/prim/scal/fun/value_of.hpp>
#include <stan/math/prim/scal/err/check_less.hpp>
#include <stan/math/prim/scal/err/check_finite.hpp>
#include <stan/math/prim/scal/err/check_nonzero_size.hpp>
#include <stan/math/prim/mat/err/check_ordered.hpp>
#include <stan/math/prim/scal/meta/return_type.hpp>
#include <stan/math/prim/arr/functor/coupled_ode_system_cvode.hpp>
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
    integrate_ode_cvode(const F& f,
                        const std::vector<T1> y0,
                        const double t0,
                        const std::vector<double>& ts,
                        const std::vector<T2>& theta,
                        const std::vector<double>& x,
                        const std::vector<int>& x_int,
                        double rel_tol,
                        double abs_tol,
                        long int max_num_steps,
                        std::ostream* msgs) {

      stan::math::check_finite("integrate_ode", "initial state", y0);
      stan::math::check_finite("integrate_ode", "initial time", t0);
      stan::math::check_finite("integrate_ode", "times", ts);
      stan::math::check_finite("integrate_ode", "parameter vector", theta);
      stan::math::check_finite("integrate_ode", "continuous data", x);

      stan::math::check_nonzero_size("integrate_ode", "times", ts);
      stan::math::check_nonzero_size("integrate_ode", "initial state", y0);
      stan::math::check_ordered("integrate_ode", "times", ts);
      stan::math::check_less("integrate_ode", "initial time", t0, ts[0]);

      // creates basic or coupled system by template specializations
      coupled_ode_system_cvode<F, T1, T2>
        coupled_system(f, y0, theta, x, x_int,
                       rel_tol, abs_tol, max_num_steps,
                       msgs);

      // first time in the vector must be time of initial state
      std::vector<double> ts_vec(ts.size() + 1);
      ts_vec[0] = t0;
      for (size_t n = 0; n < ts.size(); n++)
        ts_vec[n + 1] = ts[n];

      std::vector<std::vector<double> > y_coupled(ts_vec.size());

      // the coupled system creates the coupled initial state
      std::vector<double> initial_coupled_state
        = coupled_system.initial_state();

      coupled_system.integrate_times(ts, y_coupled);

      // remove the first state corresponding to the initial value
      y_coupled.erase(y_coupled.begin());

      // the coupled system also encapsulates the decoupling operation
      return coupled_system.decouple_states(y_coupled);
    }

  }

}

#endif
