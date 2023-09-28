#ifndef STAN_MATH_FWD_FUNCTOR_INTEGRATE_1D_HPP
#define STAN_MATH_FWD_FUNCTOR_INTEGRATE_1D_HPP

#include <stan/math/fwd/meta.hpp>
#include <stan/math/prim/functor/integrate_1d.hpp>
#include <stan/math/prim/fun/value_of.hpp>
#include <stan/math/prim/meta/forward_as.hpp>
#include <stan/math/prim/functor/apply.hpp>
#include <stan/math/fwd/functor/finite_diff.hpp>

namespace stan {
namespace math {
/**
 * Return the integral of f from a to b to the given relative tolerance
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
 * @param[in, out] msgs the print stream for warning messages
 * @param args additional arguments to pass to f
 * @return numeric integral of function f
 */
template <typename F, typename T_a, typename T_b, typename... Args,
          require_any_st_fvar<T_a, T_b, Args...> * = nullptr>
inline return_type_t<T_a, T_b, Args...> integrate_1d_impl(
    const F &f, const T_a &a, const T_b &b, double relative_tolerance,
    std::ostream *msgs, const Args &... args) {
  using FvarT = scalar_type_t<return_type_t<T_a, T_b, Args...>>;

  // Wrap integrate_1d call in a functor where the input arguments are only
  // for which tangents are needed
  auto a_val = value_of(a);
  auto b_val = value_of(b);
  auto func
      = [f, msgs, relative_tolerance, a_val, b_val](const auto &... args_var) {
          return integrate_1d_impl(f, a_val, b_val, relative_tolerance, msgs,
                                   args_var...);
        };
  FvarT ret = finite_diff(func, args...);

  // Calculate tangents w.r.t. integration bounds if needed
  if (is_fvar<T_a>::value || is_fvar<T_b>::value) {
    auto val_args = std::make_tuple(value_of(args)...);
    if (is_fvar<T_a>::value) {
      ret.d_ += math::forward_as<FvarT>(a).d_
                * math::apply(
                    [&](auto &&... tuple_args) {
                      return -f(a_val, 0.0, msgs, tuple_args...);
                    },
                    val_args);
    }
    if (is_fvar<T_b>::value) {
      ret.d_ += math::forward_as<FvarT>(b).d_
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
 * Compute the integral of the single variable function f from a to b to within
 * a specified relative tolerance. a and b can be finite or infinite.
 *
 * @tparam T_a type of first limit
 * @tparam T_b type of second limit
 * @tparam T_theta type of parameters
 * @tparam T Type of f
 *
 * @param f the functor to integrate
 * @param a lower limit of integration
 * @param b upper limit of integration
 * @param theta additional parameters to be passed to f
 * @param x_r additional data to be passed to f
 * @param x_i additional integer data to be passed to f
 * @param[in, out] msgs the print stream for warning messages
 * @param relative_tolerance relative tolerance passed to Boost quadrature
 * @return numeric integral of function f
 */
template <typename F, typename T_a, typename T_b, typename T_theta,
          require_any_fvar_t<T_a, T_b, T_theta> * = nullptr>
inline return_type_t<T_a, T_b, T_theta> integrate_1d(
    const F &f, const T_a &a, const T_b &b, const std::vector<T_theta> &theta,
    const std::vector<double> &x_r, const std::vector<int> &x_i,
    std::ostream *msgs, const double relative_tolerance) {
  return integrate_1d_impl(integrate_1d_adapter<F>(f), a, b, relative_tolerance,
                           msgs, theta, x_r, x_i);
}

}  // namespace math
}  // namespace stan
#endif
