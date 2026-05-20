#ifndef STAN_MATH_FWD_FUNCTOR_INTEGRATE_1D_DOUBLE_EXPONENTIAL_HPP
#define STAN_MATH_FWD_FUNCTOR_INTEGRATE_1D_DOUBLE_EXPONENTIAL_HPP

#include <stan/math/fwd/meta.hpp>
#include <stan/math/prim/functor/integrate_1d_double_exponential.hpp>
#include <stan/math/prim/fun/value_of.hpp>
#include <stan/math/prim/meta/forward_as.hpp>
#include <stan/math/prim/functor/apply.hpp>
#include <stan/math/fwd/functor/finite_diff.hpp>

namespace stan {
namespace math {
/**
 * Return the integral of f from a to b using adaptive double-exponential
 * quadrature, with tangents computed via finite differences over the
 * integrand parameters.
 *
 * @tparam F Type of f
 * @tparam T_a type of first limit
 * @tparam T_b type of second limit
 * @tparam Args types of parameter pack arguments
 *
 * @param f the functor to integrate
 * @param a lower limit of integration
 * @param b upper limit of integration
 * @param relative_tolerance relative tolerance passed to Boost quadrature
 * @param absolute_tolerance absolute-error floor on the convergence test
 * @param max_refinements maximum refinement level passed to the Boost
 *   quadrature class constructor
 * @param[in, out] msgs the print stream for warning messages
 * @param args additional arguments to pass to f
 * @return numeric integral of function f
 */
template <typename F, typename T_a, typename T_b, typename... Args,
          require_any_st_fvar<T_a, T_b, Args...> * = nullptr>
inline return_type_t<T_a, T_b, Args...> integrate_1d_double_exponential_impl(
    const F &f, const T_a &a, const T_b &b, double relative_tolerance,
    double absolute_tolerance, int max_refinements, std::ostream *msgs,
    const Args &... args) {
  using FvarT = scalar_type_t<return_type_t<T_a, T_b, Args...>>;

  auto a_val = value_of(a);
  auto b_val = value_of(b);
  auto func = [f, msgs, relative_tolerance, absolute_tolerance, max_refinements,
               a_val, b_val](const auto &... args_var) {
    return integrate_1d_double_exponential_impl(
        f, a_val, b_val, relative_tolerance, absolute_tolerance,
        max_refinements, msgs, args_var...);
  };
  FvarT ret = finite_diff(func, args...);

  if constexpr (is_fvar<T_a>::value || is_fvar<T_b>::value) {
    auto val_args = std::make_tuple(value_of(args)...);
    if constexpr (is_fvar<T_a>::value) {
      ret.d_ += a.d_
                * math::apply(
                    [&](auto &&... tuple_args) {
                      return -f(a_val, 0.0, msgs, tuple_args...);
                    },
                    val_args);
    }
    if constexpr (is_fvar<T_b>::value) {
      ret.d_ += b.d_
                * math::apply(
                    [&](auto &&... tuple_args) {
                      return f(b_val, 0.0, msgs, tuple_args...);
                    },
                    val_args);
    }
  }
  return ret;
}

/**
 * Compute the integral of the single variable function f from a to b using
 * adaptive double-exponential quadrature. a and b can be finite or
 * infinite.
 *
 * @tparam T_a type of first limit
 * @tparam T_b type of second limit
 * @tparam T_theta type of parameters
 * @tparam F Type of f
 *
 * @param f the functor to integrate
 * @param a lower limit of integration
 * @param b upper limit of integration
 * @param theta additional parameters to be passed to f
 * @param x_r additional data to be passed to f
 * @param x_i additional integer data to be passed to f
 * @param[in, out] msgs the print stream for warning messages
 * @param relative_tolerance relative tolerance passed to Boost quadrature
 * @param absolute_tolerance absolute-error floor on the convergence test
 * @param max_refinements maximum refinement level passed to the Boost
 *   quadrature class constructor
 * @return numeric integral of function f
 */
template <typename F, typename T_a, typename T_b, typename T_theta,
          require_any_fvar_t<T_a, T_b, T_theta> * = nullptr>
inline return_type_t<T_a, T_b, T_theta> integrate_1d_double_exponential(
    const F &f, const T_a &a, const T_b &b, const std::vector<T_theta> &theta,
    const std::vector<double> &x_r, const std::vector<int> &x_i,
    std::ostream *msgs, const double relative_tolerance,
    const double absolute_tolerance = 0.0,
    const int max_refinements
    = INTEGRATE_1D_DOUBLE_EXPONENTIAL_MAX_REFINEMENTS) {
  return integrate_1d_double_exponential_impl(
      integrate_1d_adapter<F>(f), a, b, relative_tolerance, absolute_tolerance,
      max_refinements, msgs, theta, x_r, x_i);
}

}  // namespace math
}  // namespace stan
#endif
