#ifndef STAN_MATH_REV_FUNCTOR_DAE_HPP
#define STAN_MATH_REV_FUNCTOR_DAE_HPP

#include <stan/math/prim/functor/apply.hpp>
#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/functor/idas_integrator.hpp>
#include <stan/math/rev/functor/dae_system.hpp>
#include <ostream>
#include <vector>

namespace stan {
namespace math {

/**
 * Solve the DAE initial value problem f(t, y, y')=0, y(t0) = yy0, y'(t0)=yp0 at
 * a set of times, { t1, t2, t3, ... } using IDAS.
 *
 * \p f must define an operator() with the signature as:
 *   template<typename T_yy, typename T_yp, typename... T_Args>
 *   Eigen::Matrix<stan::return_type_t<T_yy, T_yp, T_Args...>, Eigen::Dynamic,
 * 1> operator()(double t, const Eigen::Matrix<T_yy, Eigen::Dynamic, 1>& yy,
 *     const Eigen::Matrix<T_yp, Eigen::Dynamic, 1>& yp,
 *     std::ostream* msgs, const T_Args&... args);
 *
 * t is the time, yy the vector-valued state, yp the vector-valued
 * state derivative, msgs a stream for error
 * messages, and args the optional arguments passed to the DAE solve function
 * (which are passed through to \p f without modification).
 *
 * @tparam F Type of DAE residual functor
 * @tparam T_yy0 Type of initial state
 * @tparam T_yp0 Type of initial state derivatives
 * @tparam T_Args Types of pass-through parameters
 *
 * @param func Calling function name (for printing debugging messages)
 * @param f DAE residual functor
 * @param yy0 Initial state
 * @param yp0 Initial state derivatives
 * @param t0 Initial time
 * @param ts Times at which to solve the DAE at. All values must be sorted and
 *   not less than t0.
 * @param rtol Relative tolerance passed to IDAS
 * @param atol Absolute tolerance passed to IDAS
 * @param max_num_steps Upper limit on the number of integration steps to
 *   take between each output (error if exceeded)
 * @param[in, out] msgs the print stream for warning messages
 * @param args Extra arguments passed unmodified through to DAE right hand side
 * @return Solution to DAE at times \p ts
 */
template <typename F, typename T_yy, typename T_yp, typename... T_Args,
          require_all_eigen_col_vector_t<T_yy, T_yp>* = nullptr>
std::vector<Eigen::Matrix<stan::return_type_t<T_yy, T_yp, T_Args...>, -1, 1>>
dae_tol_impl(const char* func, const F& f, const T_yy& yy0, const T_yp& yp0,
             double t0, const std::vector<double>& ts, double rtol, double atol,
             int64_t max_num_steps, std::ostream* msgs, const T_Args&... args) {
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
  math::apply(
      [&](auto&&... args) {
        std::vector<int> unused_temp{
            0, (check_finite("dae", "DAE parameters and data", args), 0)...};
      },
      args_ref_tuple);

  return math::apply(
      [&](const auto&... args_refs) {
        dae_system<F, T_yy, T_yp, ref_type_t<T_Args>...> dae(f, yy0, yp0, msgs,
                                                             args_refs...);
        idas_integrator integ(rtol, atol, max_num_steps);
        return integ(func, dae, t0, ts);
      },
      args_ref_tuple);
}

/**
 * Solve the DAE initial value problem f(t, y, y')=0, y(t0) = yy0, y'(t0)=yp0 at
 * a set of times, { t1, t2, t3, ... } using IDAS.
 *
 * \p f must define an operator() with the signature as:
 *   template<typename T_yy, typename T_yp, typename... T_Args>
 *   Eigen::Matrix<stan::return_type_t<T_yy, T_yp, T_Args...>, Eigen::Dynamic,
 * 1> operator()(double t, const Eigen::Matrix<T_yy, Eigen::Dynamic, 1>& yy,
 *     const Eigen::Matrix<T_yp, Eigen::Dynamic, 1>& yp,
 *     std::ostream* msgs, const T_Args&... args);
 *
 * t is the time, yy the vector-valued state, yp the vector-valued
 * state derivative, msgs a stream for error
 * messages, and args the optional arguments passed to the DAE solve function
 * (which are passed through to \p f without modification).
 *
 * @tparam F Type of DAE residual functor
 * @tparam T_yy0 Type of initial state
 * @tparam T_yp0 Type of initial state derivatives
 * @tparam T_Args Types of pass-through parameters
 *
 * @param f DAE residual functor
 * @param yy0 Initial state
 * @param yp0 Initial state derivatives
 * @param t0 Initial time
 * @param ts Times at which to solve the DAE at. All values must be sorted and
 *   not less than t0.
 * @param rtol Relative tolerance passed to IDAS
 * @param atol Absolute tolerance passed to IDAS
 * @param max_num_steps Upper limit on the number of integration steps to
 *   take between each output (error if exceeded)
 * @param[in, out] msgs the print stream for warning messages
 * @param args Extra arguments passed unmodified through to DAE right hand side
 * @return Solution to DAE at times \p ts
 */
template <typename F, typename T_yy, typename T_yp, typename... T_Args,
          require_all_eigen_col_vector_t<T_yy, T_yp>* = nullptr>
std::vector<Eigen::Matrix<stan::return_type_t<T_yy, T_yp, T_Args...>, -1, 1>>
dae_tol(const F& f, const T_yy& yy0, const T_yp& yp0, double t0,
        const std::vector<double>& ts, double rtol, double atol,
        int64_t max_num_steps, std::ostream* msgs, const T_Args&... args) {
  return dae_tol_impl("dae_tol", f, yy0, yp0, t0, ts, rtol, atol, max_num_steps,
                      msgs, args...);
}

/**
 * Solve the DAE initial value problem f(t, y, y')=0, y(t0) = yy0, y'(t0)=yp0 at
 * a set of times, { t1, t2, t3, ... } using IDAS, assuming default controls
 * (relative tol, absolute tol, max number of steps) = (1.e-10, 1.e-10, 1e8).
 *
 * \p f must define an operator() with the signature as:
 *   template<typename T_yy, typename T_yp, typename... T_Args>
 *   Eigen::Matrix<stan::return_type_t<T_yy, T_yp, T_Args...>, Eigen::Dynamic,
 * 1> operator()(double t, const Eigen::Matrix<T_yy, Eigen::Dynamic, 1>& yy,
 *     const Eigen::Matrix<T_yp, Eigen::Dynamic, 1>& yp,
 *     std::ostream* msgs, const T_Args&... args);
 *
 * t is the time, yy the vector-valued state, yp the vector-valued
 * state derivative, msgs a stream for error
 * messages, and args the optional arguments passed to the DAE solve function
 * (which are passed through to \p f without modification).
 *
 * @tparam F Type of DAE residual functor
 * @tparam T_yy0 Type of initial state
 * @tparam T_yp0 Type of initial state derivatives
 * @tparam T_Args Types of pass-through parameters
 *
 * @param f DAE residual functor
 * @param yy0 Initial state
 * @param yp0 Initial state derivatives
 * @param t0 Initial time
 * @param ts Times at which to solve the DAE at. All values must be sorted and
 *   not less than t0.
 * @param[in, out] msgs the print stream for warning messages
 * @param args Extra arguments passed unmodified through to DAE right hand side
 * @return Solution to DAE at times \p ts
 */
template <typename F, typename T_yy, typename T_yp, typename... T_Args,
          require_all_eigen_col_vector_t<T_yy, T_yp>* = nullptr>
std::vector<Eigen::Matrix<stan::return_type_t<T_yy, T_yp, T_Args...>, -1, 1>>
dae(const F& f, const T_yy& yy0, const T_yp& yp0, double t0,
    const std::vector<double>& ts, std::ostream* msgs, const T_Args&... args) {
  return dae_tol_impl("dae", f, yy0, yp0, t0, ts, 1.e-10, 1.e-10, 1e8, msgs,
                      args...);
}

}  // namespace math
}  // namespace stan

#endif
