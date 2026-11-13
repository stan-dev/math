#ifndef STAN_MATH_PRIM_FUN_LOG_BESSEL_I0E_HPP
#define STAN_MATH_PRIM_FUN_LOG_BESSEL_I0E_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/bessel_i0e_derivative_ratio.hpp>
#include <stan/math/prim/fun/log_modified_bessel_first_kind.hpp>
#include <stan/math/prim/fun/value_of.hpp>
#include <stan/math/prim/functor/apply_scalar_unary.hpp>
#include <stan/math/prim/functor/partials_propagator.hpp>
#include <boost/math/tools/rational.hpp>
#include <cmath>
#include <utility>

namespace stan {
namespace math {

/**
 * Return the log of the exponentially scaled modified Bessel function of the
 * first kind of order zero, \f$\log(I_0(z) e^{-z}) = \log I_0(z) - z\f$.
 *
 * For \f$z \ge 7.75\f$ the value comes from Boost's approximation of
 * \f$\sqrt{z} e^{-z} I_0(z)\f$ (the coefficients of
 * `log_modified_bessel_first_kind`) without adding and subtracting \f$z\f$,
 * so it keeps its absolute precision for every \f$z\f$. The derivative is
 * \f$-(1 - I_1(z) / I_0(z))\f$, from `bessel_i0e_derivative_ratio`, which
 * also supplies the derivatives of higher order.
 *
 * @tparam T type of the argument
 * @param z non-negative argument
 * @return \f$\log I_0(z) - z\f$
 * @throw std::domain_error if `z` is NaN or negative
 */
template <typename T, require_stan_scalar_t<T>* = nullptr,
          require_not_complex_t<T>* = nullptr>
inline return_type_t<T> log_bessel_i0e(T&& z) {
  if constexpr (is_autodiff_v<T>) {
    const auto z_val = value_of(z);
    const auto log_i0e = log_bessel_i0e(z_val);
    auto ops_partials = make_partials_propagator(z);
    partials<0>(ops_partials)[0] = -bessel_i0e_derivative_ratio(0, z_val);
    return ops_partials.build(log_i0e);
  } else {
    const double x = z;
    check_nonnegative("log_bessel_i0e", "z", x);
    using boost::math::tools::evaluate_polynomial;
    if (x < 7.75) {
      const double q = 0.25 * x * x;
      return std::log1p(q * evaluate_polynomial(internal::bessel_i0_small_p, q))
             - x;
    }
    if (x < 500) {
      return std::log(
                 evaluate_polynomial(internal::bessel_i0_medium_p, 1.0 / x))
             - 0.5 * std::log(x);
    }
    return std::log(evaluate_polynomial(internal::bessel_i0_large_p, 1.0 / x))
           - 0.5 * std::log(x);
  }
}

/**
 * Structure to wrap `log_bessel_i0e()` so it can be vectorized.
 */
struct log_bessel_i0e_fun {
  /**
   * Return the log of the exponentially scaled Bessel function \f$I_0\f$.
   *
   * @tparam T type of the argument
   * @param x argument
   * @return \f$\log I_0(x) - x\f$
   */
  template <typename T>
  static inline auto fun(T&& x) {
    return log_bessel_i0e(std::forward<T>(x));
  }
};

/**
 * Return the elementwise `log_bessel_i0e()` of a container.
 *
 * @tparam T type of the container
 * @param x container
 * @return \f$\log I_0(x) - x\f$ of each element
 */
template <
    typename T,
    require_all_not_nonscalar_prim_or_rev_kernel_expression_t<T>* = nullptr,
    require_container_t<T>* = nullptr, require_not_var_matrix_t<T>* = nullptr>
inline auto log_bessel_i0e(T&& x) {
  return apply_scalar_unary<log_bessel_i0e_fun, T>::apply(std::forward<T>(x));
}

}  // namespace math
}  // namespace stan

#endif
