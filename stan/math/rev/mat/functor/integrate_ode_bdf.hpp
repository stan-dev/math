#ifndef STAN_MATH_REV_MAT_FUNCTOR_INTEGRATE_ODE_BDF_HPP
#define STAN_MATH_REV_MAT_FUNCTOR_INTEGRATE_ODE_BDF_HPP

#include <stan/math/prim/arr/fun/value_of.hpp>
#include <stan/math/prim/scal/err/check_less.hpp>
#include <stan/math/prim/scal/err/check_finite.hpp>
#include <stan/math/prim/arr/err/check_nonzero_size.hpp>
#include <stan/math/prim/arr/err/check_ordered.hpp>
#include <stan/math/prim/scal/meta/return_type.hpp>
#include <stan/math/rev/mat/functor/cvodes_integrator.hpp>
#include <stan/math/rev/mat/functor/ode_system.hpp>
#include <stan/math/rev/arr/fun/decouple_ode_states.hpp>
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
     * This function is templated to allow the initial times to be
     * either data or autodiff variables and the parameters to be data
     * or autodiff variables.  The autodiff-based implementation for
     * reverse-mode are defined in namespace <code>stan::math</code>
     * and may be invoked via argument-dependent lookup by including
     * their headers.
     *
     * The solver used is based on the backward differentiation
     * formula which is an implicit numerical integration scheme
     * appropiate for stiff ODE systems.
     *
     * @tparam F type of ODE system function.
     * @tparam T_initial type of scalars for initial values.
     * @tparam T_param type of scalars for parameters.
     * @param[in] f functor for the base ordinary differential equation.
     * @param[in] y0 initial state.
     * @param[in] t0 initial time.
     * @param[in] ts times of the desired solutions, in strictly
     * increasing order, all greater than the initial time.
     * @param[in] theta parameter vector for the ODE.
     * @param[in] x continuous data vector for the ODE.
     * @param[in] x_int integer data vector for the ODE.
     * @param[in] rel_tol relative tolerance of solution
     * @param[in] abs_tol absolute tolerance of solution
     * @param[in] max_num_steps maximal number of admissable steps
     * between time-points
     * @param[in, out] msgs the print stream for warning messages.
     * @return a vector of states, each state being a vector of the
     * same size as the state variable, corresponding to a time in ts.
     */
    template <typename F, typename T_initial, typename T_param>
    std::vector<std::vector<typename stan::return_type<T_initial,
                                                       T_param>::type> >
    integrate_ode_bdf(const F& f,
                      const std::vector<T_initial>& y0,
                      const double t0,
                      const std::vector<double>& ts,
                      const std::vector<T_param> theta,
                      const std::vector<double>& x,
                      const std::vector<int>& x_int,
                      double rel_tol = 1e-10,
                      double abs_tol = 1e-10,
                      long int max_num_steps = 1e8,  // NOLINT(runtime/int)
                      std::ostream* msgs = 0) {
      stan::math::check_finite("integrate_ode_bdf",
                               "initial state", y0);
      stan::math::check_finite("integrate_ode_bdf",
                               "initial time", t0);
      stan::math::check_finite("integrate_ode_bdf",
                               "times", ts);
      stan::math::check_finite("integrate_ode_bdf",
                               "parameter vector", theta);
      stan::math::check_finite("integrate_ode_bdf",
                               "continuous data", x);

      stan::math::check_nonzero_size("integrate_ode_bdf",
                                     "times", ts);
      stan::math::check_nonzero_size("integrate_ode_bdf",
                                     "initial state", y0);
      stan::math::check_ordered("integrate_ode_bdf",
                                "times", ts);
      stan::math::check_less("integrate_ode_bdf",
                             "initial time", t0, ts[0]);

      cvodes_integrator<F, T_initial, T_param>
        integrator(f, y0, t0, theta,
                   x, x_int,
                   ts,
                   rel_tol, abs_tol,
                   max_num_steps,
                   1,
                   msgs);

      return integrator.integrate();
    }
  }  // math
}  // stan
#endif
