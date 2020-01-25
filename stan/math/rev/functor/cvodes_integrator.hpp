#ifndef STAN_MATH_REV_FUNCTOR_INTEGRATE_ODE_CVODES_HPP
#define STAN_MATH_REV_FUNCTOR_INTEGRATE_ODE_CVODES_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/functor/cvodes_utils.hpp>
#include <stan/math/rev/functor/coupled_ode_system.hpp>
#include <stan/math/rev/functor/cvodes_ode_data.hpp>
#include <stan/math/rev/functor/coupled_ode_observer.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/value_of.hpp>
#include <cvodes/cvodes.h>
#include <sunlinsol/sunlinsol_dense.h>
#include <algorithm>
#include <ostream>
#include <vector>

namespace stan {
namespace math {

/**
 * Integrator interface for CVODES' ODE solvers (Adams & BDF
 * methods).
 *
 * @tparam Lmm ID of ODE solver (1: ADAMS, 2: BDF)
 */
template <int Lmm>
class cvodes_integrator {
 public:
  cvodes_integrator() {}

  /**
   * Return the solutions for the specified system of ordinary
   * differential equations given the specified initial state,
   * initial times, times of desired solution, and parameters and
   * data, writing error and warning messages to the specified
   * stream.
   *
   * This function is templated to allow the initials to be
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
   * @tparam T_t0 type of scalar of initial time point.
   * @tparam T_ts type of time-points where ODE solution is returned.
   *
   * @param[in] f functor for the base ordinary differential equation.
   * @param[in] y0 initial state.
   * @param[in] t0 initial time.
   * @param[in] ts times of the desired solutions, in strictly
   * increasing order, all greater than the initial time.
   * @param[in] theta parameter vector for the ODE.
   * @param[in] x continuous data vector for the ODE.
   * @param[in] x_int integer data vector for the ODE.
   * @param[in, out] msgs the print stream for warning messages.
   * @param[in] relative_tolerance relative tolerance passed to CVODE.
   * @param[in] absolute_tolerance absolute tolerance passed to CVODE.
   * @param[in] max_num_steps maximal number of admissable steps
   * between time-points
   * @return a vector of states, each state being a vector of the
   * same size as the state variable, corresponding to a time in ts.
   */
  template <typename... Args, typename F, typename T_initial,
	    typename T_t0, typename T_ts>
  std::vector<std::vector<
		typename stan::return_type<T_initial, T_t0, T_ts, Args...>::type>>
  integrate(const F& f, const std::vector<T_initial>& y0, const T_t0& t0,
            const std::vector<T_ts>& ts, const Args&... args,
            std::ostream* msgs, double relative_tolerance,
            double absolute_tolerance,
            long int max_num_steps) {  // NOLINT(runtime/int)
    using initial_var = stan::is_var<T_initial>;

    const char* fun = "integrate_ode_cvodes";

    const double t0_dbl = value_of(t0);
    const std::vector<double> ts_dbl = value_of(ts);

    check_finite(fun, "initial state", y0);
    check_finite(fun, "initial time", t0_dbl);
    check_finite(fun, "times", ts_dbl);

    // Code from: https://stackoverflow.com/a/17340003 . Should probably do something better
    std::vector<int> unused_temp{ 0, (check_finite(fun, "ode parameters and data", args), 0)... };

    check_nonzero_size(fun, "times", ts);
    check_nonzero_size(fun, "initial state", y0);
    check_ordered(fun, "times", ts_dbl);
    check_less(fun, "initial time", t0_dbl, ts_dbl[0]);
    if (relative_tolerance <= 0) {
      invalid_argument("integrate_ode_cvodes", "relative_tolerance,",
                       relative_tolerance, "", ", must be greater than 0");
    }
    if (absolute_tolerance <= 0) {
      invalid_argument("integrate_ode_cvodes", "absolute_tolerance,",
                       absolute_tolerance, "", ", must be greater than 0");
    }
    if (max_num_steps <= 0) {
      invalid_argument("integrate_ode_cvodes", "max_num_steps,", max_num_steps,
                       "", ", must be greater than 0");
    }


    const size_t y0_vars = internal::count_vars(y0);
    const size_t args_vars = internal::count_vars(args...);
    const size_t S = y0_vars + args_vars;

    using ode_data = cvodes_ode_data<F, T_initial, Args...>;
    ode_data cvodes_data(f, y0, args..., msgs);

    void* cvodes_mem = CVodeCreate(Lmm);
    if (cvodes_mem == nullptr) {
      throw std::runtime_error("CVodeCreate failed to allocate memory");
    }

    const size_t coupled_size = cvodes_data.coupled_ode_.size();

    using return_t = return_type_t<T_initial, T_t0, T_ts, Args...>;
    
    std::vector<std::vector<return_t>> y;
    coupled_ode_observer<T_initial, T_t0, T_ts, F, Args...> observer(f, y0, t0, ts, args..., msgs, y);

    try {
      check_flag_sundials(CVodeInit(cvodes_mem, &ode_data::cv_rhs, t0_dbl,
                                    cvodes_data.nv_state_),
                          "CVodeInit");

      // Assign pointer to this as user data
      check_flag_sundials(
          CVodeSetUserData(cvodes_mem, reinterpret_cast<void*>(&cvodes_data)),
          "CVodeSetUserData");

      cvodes_set_options(cvodes_mem, relative_tolerance, absolute_tolerance,
                         max_num_steps);

      // for the stiff solvers we need to reserve additional memory
      // and provide a Jacobian function call. new API since 3.0.0:
      // create matrix object and linear solver object; resource
      // (de-)allocation is handled in the cvodes_ode_data
      check_flag_sundials(
          CVodeSetLinearSolver(cvodes_mem, cvodes_data.LS_, cvodes_data.A_),
          "CVodeSetLinearSolver");
      check_flag_sundials(
          CVodeSetJacFn(cvodes_mem, &ode_data::cv_jacobian_states),
          "CVodeSetJacFn");

      // initialize forward sensitivity system of CVODES as needed
      if (S > 0) {
        check_flag_sundials(
            CVodeSensInit(cvodes_mem, static_cast<int>(S), CV_STAGGERED,
                          &ode_data::cv_rhs_sens, cvodes_data.nv_state_sens_),
            "CVodeSensInit");

        check_flag_sundials(CVodeSensEEtolerances(cvodes_mem),
                            "CVodeSensEEtolerances");
      }

      double t_init = t0_dbl;
      for (size_t n = 0; n < ts.size(); ++n) {
        double t_final = ts_dbl[n];
        if (t_final != t_init) {
          check_flag_sundials(CVode(cvodes_mem, t_final, cvodes_data.nv_state_,
                                    &t_init, CV_NORMAL),
                              "CVode");
        }
        if (S > 0) {
          check_flag_sundials(
              CVodeGetSens(cvodes_mem, &t_init, cvodes_data.nv_state_sens_),
              "CVodeGetSens");
        }
        observer(cvodes_data.coupled_state_, t_final);
        t_init = t_final;
      }
    } catch (const std::exception& e) {
      CVodeFree(&cvodes_mem);
      throw;
    }

    CVodeFree(&cvodes_mem);

    return y;
  }
};  // cvodes integrator

}  // namespace math
}  // namespace stan
#endif
