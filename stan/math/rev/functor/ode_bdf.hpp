#ifndef STAN_MATH_REV_FUNCTOR_ODE_BDF_HPP
#define STAN_MATH_REV_FUNCTOR_ODE_BDF_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/functor/cvodes_integrator.hpp>
#include <ostream>
#include <vector>

namespace stan {
namespace math {

template <typename F, typename T_initial, typename T_t0, typename T_ts,
          typename... T_Args>
std::vector<Eigen::Matrix<stan::return_type_t<T_initial, T_t0, T_ts, T_Args...>,
                          Eigen::Dynamic, 1>>
ode_bdf(const F& f, const Eigen::Matrix<T_initial, Eigen::Dynamic, 1>& y0,
        const T_t0& t0, const std::vector<T_ts>& ts, std::ostream* msgs,
        const T_Args&... args) {
  double relative_tolerance = 1e-6;
  double absolute_tolerance = 1e-6;
  long int max_num_steps = 1e6;

  return ode_bdf_tol(f, y0, t0, ts, relative_tolerance, absolute_tolerance,
                     max_num_steps, msgs, args...);
}

template <typename F, typename T_initial, typename T_t0, typename T_ts,
          typename... T_Args>
std::vector<Eigen::Matrix<stan::return_type_t<T_initial, T_t0, T_ts, T_Args...>,
                          Eigen::Dynamic, 1>>
ode_bdf_tol(const F& f, const Eigen::Matrix<T_initial, Eigen::Dynamic, 1>& y0,
            const T_t0& t0, const std::vector<T_ts>& ts,
            double relative_tolerance, double absolute_tolerance,
            long int max_num_steps, std::ostream* msgs, const T_Args&... args) {
  stan::math::cvodes_integrator<CV_BDF, F, T_initial, T_t0, T_ts, T_Args...>
  integrator(f, y0, t0, ts, relative_tolerance, absolute_tolerance,
             max_num_steps, msgs, args...);

  return integrator.integrate();
}

}  // namespace math
}  // namespace stan
#endif
