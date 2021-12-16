#ifndef STAN_MATH_REV_FUNCTOR_INTEGRATOR_DAE_HPP
#define STAN_MATH_REV_FUNCTOR_INTEGRATOR_DAE_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/functor/idas_integrator.hpp>
#include <stan/math/rev/functor/idas_system.hpp>
#include <ostream>
#include <vector>

namespace stan {
namespace math {

template <typename F, typename T_yy, typename T_yp,
          typename... T_Args, require_all_eigen_col_vector_t<T_yy, T_yp>* = nullptr>
std::vector<Eigen::Matrix<stan::return_type_t<T_yy, T_yp, T_Args...>, -1, 1>>
dae_tol_impl(const char* func, const F& f,
             const T_yy& yy0, const T_yp& yp0,
             double t0, const std::vector<double>& ts,
             double rtol, double atol, long int max_num_steps,
             std::ostream* msgs, const T_Args&... args) {
  check_finite(func, "initial state", yy0);
  check_finite(func, "initial state derivative", yp0);
  check_nonzero_size(func, "initial state", yy0);
  check_nonzero_size(func, "initial state derivative", yp0);
  check_finite(func, "initial time", t0);
  check_finite(func, "times", ts);
  check_nonzero_size(func, "times", ts);
  check_sorted(func, "times", ts);
  check_less(func, "initial time", t0, ts[0]);
  check_positive_finite(func, "relative_tolerance", rtol);
  check_positive_finite(func, "absolute_tolerance", atol);
  check_positive(func, "max_num_steps", max_num_steps);

  const auto& args_ref_tuple = std::make_tuple(to_ref(args)...);
  return apply([&](const auto&... args_refs) {
                 dae_system<F, T_yy, T_yp, ref_type_t<T_Args>...> dae(f, yy0, yp0, msgs, args_refs...);
                 idas_integrator integ(rtol, atol, max_num_steps);
                 return integ(func, dae, t0, ts);
               },
               args_ref_tuple);
}

/**
 * Solve the ODE initial value problem y' = f(t, y), y(t0) = y0 at a set of
 * times, { t1, t2, t3, ... } using the stiff backward differentiation formula
 * BDF solver from CVODES.
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
template <typename F, typename T_yy, typename T_yp,
          typename... T_Args, require_all_eigen_col_vector_t<T_yy, T_yp>* = nullptr>
std::vector<Eigen::Matrix<stan::return_type_t<T_yy, T_yp, T_Args...>, -1, 1>>
dae_tol(const F& f, const T_yy& yy0, const T_yp& yp0,
        double t0, const std::vector<double>& ts,
        double rtol, double atol, long int max_num_steps,
        std::ostream* msgs, const T_Args&... args) {
  return dae_tol_impl("dae_tol", f, yy0, yp0, t0, ts, rtol,
                      atol, max_num_steps, msgs, args...);
}

/**
 * Return the solutions for a semi-explicit DAE system with residual
 * specified by functor F,
 * given the specified consistent initial state yy0 and yp0.
 *
 * @tparam DAE type of DAE system
 * @tparam Tpar scalar type of parameter theta
 *
 * @param[in] f functor for the base ordinary differential equation
 * @param[in] yy0 initial state
 * @param[in] yp0 initial derivative state
 * @param[in] t0 initial time
 * @param[in] ts times of the desired solutions, in strictly
 * increasing order, all greater than the initial time
 * @param[in] theta parameters
 * @param[in] x_r real data
 * @param[in] x_i int data
 * @param[in] rtol relative tolerance passed to IDAS, required <10^-3
 * @param[in] atol absolute tolerance passed to IDAS, problem-dependent
 * @param[in] max_num_steps maximal number of admissable steps
 * between time-points
 * @param[in] msgs message
 * @return a vector of states, each state being a vector of the
 * same size as the state variable, corresponding to a time in ts.
 */
template <typename F, typename T_yy, typename T_yp,
          typename... T_Args, require_all_eigen_col_vector_t<T_yy, T_yp>* = nullptr>
std::vector<Eigen::Matrix<stan::return_type_t<T_yy, T_yp, T_Args...>, -1, 1>>
dae(const F& f, const T_yy& yy0, const T_yp& yp0,
    double t0, const std::vector<double>& ts,
    double rtol, double atol, long int max_num_steps,
    std::ostream* msgs, const T_Args&... args) {
  return dae_tol_impl("dae", f, yy0, yp0, t0, ts, 
                      1.e-10, 1.e-10, 1e8, msgs, args...);
}

}  // namespace math
}  // namespace stan

#endif
