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
std::vector<std::vector<
    typename stan::return_type<T_initial, T_t0, T_ts, T_Args...>::type>>
ode_bdf(const F& f, const std::vector<T_initial>& y0, const T_t0& t0,
        const std::vector<T_ts>& ts, double relative_tolerance,
        double absolute_tolerance, long int max_num_steps, std::ostream* msgs,
        const T_Args&... args) {  // NOLINT(runtime/int)
  stan::math::cvodes_integrator<CV_BDF, F, T_initial, T_t0, T_ts, T_Args...>
  integrator(f, y0, t0, ts, args..., msgs, relative_tolerance,
             absolute_tolerance, max_num_steps);
  return integrator.integrate();
}

}  // namespace math
}  // namespace stan
#endif
