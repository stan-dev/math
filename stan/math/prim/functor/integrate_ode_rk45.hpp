#ifndef STAN_MATH_PRIM_FUNCTOR_INTEGRATE_ODE_RK45_HPP
#define STAN_MATH_PRIM_FUNCTOR_INTEGRATE_ODE_RK45_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/functor/closure_adapter.hpp>
#include <stan/math/prim/functor/integrate_ode_std_vector_interface_adapter.hpp>
#include <stan/math/prim/functor/ode_rk45.hpp>
#include <ostream>
#include <vector>

namespace stan {
namespace math {

namespace internal {

template <typename F, typename T_y0, typename T_param, typename T_t0,
          typename T_ts, require_not_stan_closure_t<F>* = nullptr>
inline auto integrate_ode_rk45_impl(
    const F& f, const std::vector<T_y0>& y0, const T_t0& t0,
    const std::vector<T_ts>& ts, const std::vector<T_param>& theta,
    const std::vector<double>& x, const std::vector<int>& x_int,
    std::ostream* msgs, double relative_tolerance, double absolute_tolerance,
    int max_num_steps) {
  internal::integrate_ode_std_vector_interface_adapter<F> f_adapted(f);
  return ode_rk45_tol_impl("integrate_ode_rk45", f_adapted, to_vector(y0), t0,
                           ts, relative_tolerance, absolute_tolerance,
                           max_num_steps, msgs, theta, x, x_int);
}

template <typename F, typename T_y0, typename T_param, typename T_t0,
          typename T_ts, require_stan_closure_t<F>* = nullptr>
inline auto integrate_ode_rk45_impl(
    const F& f, const std::vector<T_y0>& y0, const T_t0& t0,
    const std::vector<T_ts>& ts, const std::vector<T_param>& theta,
    const std::vector<double>& x, const std::vector<int>& x_int,
    std::ostream* msgs, double relative_tolerance, double absolute_tolerance,
    int max_num_steps) {
  return ode_rk45_tol_impl("integrate_ode_rk45",
                           integrate_ode_closure_adapter(), to_vector(y0), t0,
                           ts, relative_tolerance, absolute_tolerance,
                           max_num_steps, msgs, f, theta, x, x_int);
}

}  // namespace internal

/**
 * @deprecated use <code>ode_rk45</code>
 */
template <typename F, typename T_y0, typename T_param, typename T_t0,
          typename T_ts>
inline auto integrate_ode_rk45(
    const F& f, const std::vector<T_y0>& y0, const T_t0& t0,
    const std::vector<T_ts>& ts, const std::vector<T_param>& theta,
    const std::vector<double>& x, const std::vector<int>& x_int,
    std::ostream* msgs = nullptr, double relative_tolerance = 1e-6,
    double absolute_tolerance = 1e-6, int max_num_steps = 1e6) {
  auto y = internal::integrate_ode_rk45_impl(f, y0, t0, ts, theta, x, x_int,
                                             msgs, relative_tolerance,
                                             absolute_tolerance, max_num_steps);

  std::vector<std::vector<fn_return_type_t<F, T_y0, T_param, T_t0, T_ts>>>
      y_converted;
  y_converted.reserve(y.size());
  for (size_t i = 0; i < y.size(); ++i)
    y_converted.emplace_back(to_array_1d(y[i]));

  return y_converted;
}

}  // namespace math
}  // namespace stan

#endif
