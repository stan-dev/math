#ifndef STAN_MATH_PRIM_FUNCTOR_ODE_RK45_HPP
#define STAN_MATH_PRIM_FUNCTOR_ODE_RK45_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/functor/apply.hpp>
#include <stan/math/prim/functor/arkode_integrator.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <ostream>
#include <tuple>
#include <vector>

namespace stan {
namespace math {

/**
 * Solve the ODE initial value problem y' = f(t, y), y(t0) = y0 at a set of
 * times, { t1, t2, t3, ... } using SUNDIALS ARKODE's ERKStep non-stiff
 * Dormand-Prince 5(4) solver.
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
 * t is the time, y is the vector-valued state, msgs is a stream for error
 * messages, and args are optional arguments passed to the ODE solve function
 * (which are passed through to \p f without modification).
 *
 * @tparam F Type of ODE right hand side
 * @tparam T_y0 Type of initial condition
 * @tparam T_t0 Type of initial time
 * @tparam T_ts Type of output times
 * @tparam T_Args Types of pass-through parameters
 *
 * @param function_name Calling function name (for printing debugging messages)
 * @param f Right hand side of the ODE
 * @param y0_arg Initial state
 * @param t0 Initial time
 * @param ts Times at which to solve the ODE at. All values must be sorted and
 *   greater than t0.
 * @param relative_tolerance Relative tolerance passed to ARKODE
 * @param absolute_tolerance Absolute tolerance passed to ARKODE
 * @param max_num_steps Upper limit on the number of integration steps to
 *   take between each output (error if exceeded)
 * @param[in, out] msgs the print stream for warning messages
 * @param args Extra arguments passed unmodified through to ODE right hand side
 * @return Solution to ODE at times \p ts
 */
template <typename F, typename T_y0, typename T_t0, typename T_ts,
          typename... Args, require_eigen_vector_t<T_y0>* = nullptr>
inline std::vector<Eigen::Matrix<stan::return_type_t<T_y0, T_t0, T_ts, Args...>,
                                 Eigen::Dynamic, 1>>
ode_rk45_tol_impl(const char* function_name, const F& f, const T_y0& y0_arg,
                  T_t0 t0, const std::vector<T_ts>& ts,
                  double relative_tolerance, double absolute_tolerance,
                  long int max_num_steps,  // NOLINT(runtime/int)
                  std::ostream* msgs, const Args&... args) {
  const auto& args_ref_tuple = std::make_tuple(to_ref(args)...);
  return math::apply(
      [&](const auto&... args_refs) {
        arkode_integrator<ARKODE_DORMAND_PRINCE_7_4_5, F, T_y0, T_t0, T_ts,
                          ref_type_t<Args>...>
            integrator(function_name, f, y0_arg, t0, ts, relative_tolerance,
                       absolute_tolerance, max_num_steps, msgs, args_refs...);

        return integrator();
      },
      args_ref_tuple);
}

/**
 * Solve the ODE initial value problem y' = f(t, y), y(t0) = y0 at a set of
 * times, { t1, t2, t3, ... } using SUNDIALS ARKODE's ERKStep non-stiff
 * Dormand-Prince 5(4) solver.
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
 * @param y0_arg Initial state
 * @param t0 Initial time
 * @param ts Times at which to solve the ODE at. All values must be sorted and
 *   greater than t0.
 * @param relative_tolerance Relative tolerance passed to ARKODE
 * @param absolute_tolerance Absolute tolerance passed to ARKODE
 * @param max_num_steps Upper limit on the number of integration steps to
 *   take between each output (error if exceeded)
 * @param[in, out] msgs the print stream for warning messages
 * @param args Extra arguments passed unmodified through to ODE right hand side
 * @return Solution to ODE at times \p ts
 */
template <typename F, typename T_y0, typename T_t0, typename T_ts,
          typename... Args, require_eigen_vector_t<T_y0>* = nullptr>
inline std::vector<Eigen::Matrix<stan::return_type_t<T_y0, T_t0, T_ts, Args...>,
                                 Eigen::Dynamic, 1>>
ode_rk45_tol(const F& f, const T_y0& y0_arg, T_t0 t0,
             const std::vector<T_ts>& ts, double relative_tolerance,
             double absolute_tolerance,
             long int max_num_steps,  // NOLINT(runtime/int)
             std::ostream* msgs, const Args&... args) {
  return ode_rk45_tol_impl("ode_rk45_tol", f, y0_arg, t0, ts,
                           relative_tolerance, absolute_tolerance,
                           max_num_steps, msgs, args...);
}

/**
 * Solve the ODE initial value problem y' = f(t, y), y(t0) = y0 at a set of
 * times, { t1, t2, t3, ... } using SUNDIALS ARKODE's ERKStep non-stiff
 * Dormand-Prince 5(4) solver with defaults for relative_tolerance,
 * absolute_tolerance, and max_num_steps.
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
 * t is the time, y is the vector-valued state, msgs is a stream for error
 * messages, and args are optional arguments passed to the ODE solve function
 * (which are passed through to \p f without modification).
 *
 * @tparam F Type of ODE right hand side
 * @tparam T_y0 Type of initial time
 * @tparam T_ts Type of output times
 * @tparam T_Args Types of pass-through parameters
 *
 * @param f Right hand side of the ODE
 * @param y0 Initial state
 * @param t0 Initial time
 * @param ts Times at which to solve the ODE at. All values must be sorted and
 *   greater than t0.
 * @param[in, out] msgs the print stream for warning messages
 * @param args Extra arguments passed unmodified through to ODE right hand side
 * @return Solution to ODE at times \p ts
 */
template <typename F, typename T_y0, typename T_t0, typename T_ts,
          typename... Args, require_eigen_vector_t<T_y0>* = nullptr>
inline std::vector<Eigen::Matrix<stan::return_type_t<T_y0, T_t0, T_ts, Args...>,
                                 Eigen::Dynamic, 1>>
ode_rk45(const F& f, const T_y0& y0, T_t0 t0, const std::vector<T_ts>& ts,
         std::ostream* msgs, const Args&... args) {
  double relative_tolerance = 1e-6;
  double absolute_tolerance = 1e-6;
  long int max_num_steps = 1e6;  // NOLINT(runtime/int)

  return ode_rk45_tol_impl("ode_rk45", f, y0, t0, ts, relative_tolerance,
                           absolute_tolerance, max_num_steps, msgs, args...);
}

}  // namespace math
}  // namespace stan
#endif
