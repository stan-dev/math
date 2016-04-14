#ifndef STAN_MATH_REV_MAT_FUNCTOR_INTEGRATE_ODE_CVODES_HPP
#define STAN_MATH_REV_MAT_FUNCTOR_INTEGRATE_ODE_CVODES_HPP

#include <stan/math/prim/scal/fun/value_of.hpp>
#include <stan/math/rev/scal/fun/value_of.hpp>
#include <stan/math/prim/scal/err/check_less.hpp>
#include <stan/math/prim/scal/err/check_bounded.hpp>
#include <stan/math/prim/scal/err/check_finite.hpp>
#include <stan/math/prim/arr/err/check_nonzero_size.hpp>
#include <stan/math/prim/arr/err/check_ordered.hpp>
#include <stan/math/prim/scal/meta/return_type.hpp>
#include <stan/math/rev/mat/functor/cvodes_integrator.hpp>
#include <stan/math/rev/arr/fun/decouple_ode_states.hpp>
#include <boost/type_traits/is_same.hpp>
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
     * The solver argument selects the solver to be used:
     * - 0 is a non-stiff Adams-Moulton
     * - 1 is a stiff BDF
     * - 2 is a stiff BDF with stability limit detection
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
     * @param[in] rel_tol relative tolerance of solution
     * @param[in] abs_tol absolute tolerance of solution
     * @param[in] max_num_steps maximal number of admissable steps
     * between time-points
     * @param[in] solver selects solver, admissable values 0,1,2
     * @param[in, out] msgs the print stream for warning messages.
     * @return a vector of states, each state being a vector of the
     * same size as the state variable, corresponding to a time in ts.
     */
    template <typename F, typename T1, typename T2>
    std::vector<std::vector<typename stan::return_type<T1, T2>::type> >
    integrate_ode_cvodes(const F& f,
                         const std::vector<T1>& y0,
                         const double t0,
                         const std::vector<double>& ts,
                         const std::vector<T2> theta,
                         const std::vector<double>& x,
                         const std::vector<int>& x_int,
                         double rel_tol = 1e-10,
                         double abs_tol = 1e-10,
                         long int max_num_steps = 1e8,  // NOLINT(runtime/int)
                         const size_t solver = 1,
                         std::ostream* msgs = 0) {
      stan::math::check_finite("integrate_ode_cvodes",
                               "initial state", y0);
      stan::math::check_finite("integrate_ode_cvodes",
                               "initial time", t0);
      stan::math::check_finite("integrate_ode_cvodes",
                               "times", ts);
      stan::math::check_finite("integrate_ode_cvodes",
                               "parameter vector", theta);
      stan::math::check_finite("integrate_ode_cvodes",
                               "continuous data", x);

      stan::math::check_nonzero_size("integrate_ode_cvodes",
                                     "times", ts);
      stan::math::check_nonzero_size("integrate_ode_cvodes",
                                     "initial state", y0);
      stan::math::check_ordered("integrate_ode_cvodes",
                                "times", ts);
      stan::math::check_less("integrate_ode_cvodes",
                             "initial time", t0, ts[0]);

      stan::math::check_bounded("integrate_ode_cvodes",
                                "solver", solver,
                                static_cast<size_t>(0),
                                static_cast<size_t>(2));

      const size_t N = y0.size();
      const size_t M = theta.size();

      typedef boost::is_same<T1, stan::math::var> initial_var;
      typedef boost::is_same<T2, stan::math::var> param_var;

      std::vector<double>    y0_dbl(N);
      std::vector<double> theta_dbl(M);

      for (size_t i = 0; i < N; i++)
        y0_dbl[i] = stan::math::value_of(y0[i]);

      for (size_t i = 0; i < M; i++)
        theta_dbl[i] = stan::math::value_of(theta[i]);

      const ode_model<F> ode(f, theta_dbl, x, x_int, msgs);

      cvodes_integrator<F> integrator(ode,
                                      y0_dbl, t0,
                                      initial_var::value,
                                      param_var::value,
                                      rel_tol, abs_tol,
                                      max_num_steps,
                                      solver);

      std::vector<std::vector<double> >
        y_res(ts.size(),
              std::vector<double>(integrator.size(), 0));

      integrator.integrate_times(ts, y_res);

      return decouple_ode_states(y_res, y0, theta);
    }
  }  // math
}  // stan
#endif
