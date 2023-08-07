#ifndef STAN_MATH_FWD_FUNCTOR_INTEGRATE_1D_HPP
#define STAN_MATH_FWD_FUNCTOR_INTEGRATE_1D_HPP

#include <stan/math/fwd/meta.hpp>
#include <stan/math/prim/functor/integrate_1d.hpp>
#include <stan/math/fwd/functor/fvar_wrapper.hpp>

namespace stan {
namespace math {

template <typename F, typename T_a, typename T_b, typename... Args,
          require_any_st_fvar<T_a, T_b, Args...> * = nullptr>
inline return_type_t<T_a, T_b, Args...> integrate_1d_impl(
    const F &f, const T_a &a, const T_b &b, double relative_tolerance,
    std::ostream *msgs, const Args &... args) {
  auto func = [&f,&msgs,&relative_tolerance,](auto&& a_var, auto&& b_var, auto&&... args_var) {
    return integrate_1d_impl(f, a_var, b_var, relative_tolerance,
                           msgs, args_var...);
  };
  return fvar_wrapper(func, a, b, args...);
}

template <typename F, typename T_a, typename T_b, typename T_theta,
          typename = require_any_fvar_t<T_a, T_b, T_theta>>
inline return_type_t<T_a, T_b, T_theta> integrate_1d(
    const F &f, const T_a &a, const T_b &b, const std::vector<T_theta> &theta,
    const std::vector<double> &x_r, const std::vector<int> &x_i,
    std::ostream *msgs, const double relative_tolerance = std::sqrt(EPSILON)) {
  auto func = [&f,&msgs,&x_r,&x_i](auto&& a_var, auto&& b_var, auto&& theta_var) {
    return integrate_1d_impl(integrate_1d_adapter<F>(f), a_var, b_var, relative_tolerance,
                           msgs, theta_var, x_r, x_i);
  };
  return fvar_wrapper(func, a, b, theta);
}

}  // namespace math
}  // namespace stan
#endif