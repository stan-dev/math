#ifndef STAN_MATH_FWD_FUN_POW_HPP
#define STAN_MATH_FWD_FUN_POW_HPP

#include <stan/math/fwd/meta.hpp>
#include <stan/math/fwd/core.hpp>
#include <stan/math/fwd/fun/inv.hpp>
#include <stan/math/fwd/fun/inv_sqrt.hpp>
#include <stan/math/fwd/fun/inv_square.hpp>
#include <stan/math/fwd/fun/log.hpp>
#include <stan/math/fwd/fun/sqrt.hpp>
#include <stan/math/fwd/fun/square.hpp>
#include <stan/math/prim/fun/pow.hpp>
#include <cmath>
#include <complex>
#include <type_traits>

namespace stan {
namespace math {
/*
 *
 * @tparam T1 Either an `fvar`, `arithmetic`, or `complex` type with an inner
 * `fvar` or `arithmetic` type.
 * @tparam T2 Either a `fvar`, `arithmetic`, or `complex` type with an inner
 * `fvar` or `arithmetic` type.
 * @param x1 Base variable.
 * @param x2 Exponent variable.
 * @return Base raised to the exponent.
 */
template <typename T1, typename T2,
          require_any_fvar_t<base_type_t<T1>, base_type_t<T2>>* = nullptr,
          require_all_stan_scalar_t<T1, T2>* = nullptr>
inline auto pow(const T1& x1, const T2& x2) {
  if constexpr (is_complex<T1>::value || is_complex<T2>::value) {
    return internal::complex_pow(x1, x2);
  } else if constexpr (is_fvar<T1>::value && is_fvar<T2>::value) {
    auto pow_x1_x2(stan::math::pow(x1.val_, x2.val_));
    return T1(pow_x1_x2,
              (x2.d_ * stan::math::log(x1.val_) + x2.val_ * x1.d_ / x1.val_)
                  * pow_x1_x2);
  } else if constexpr (is_fvar<T2>::value) {
    auto u = stan::math::pow(x1, x2.val_);
    return T2(u, x2.d_ * stan::math::log(x1) * u);
  } else {
    if (x2 == -2) {
      return stan::math::inv_square(x1);
    } else if (x2 == -1) {
      return stan::math::inv(x1);
    } else if (x2 == -0.5) {
      return stan::math::inv_sqrt(x1);
    } else if (x2 == 0.5) {
      return stan::math::sqrt(x1);
    } else if (x2 == 1.0) {
      return x1;
    } else if (x2 == 2.0) {
      return stan::math::square(x1);
    }
    return T1(stan::math::pow(x1.val_, x2),
              x1.d_ * x2 * stan::math::pow(x1.val_, x2 - 1));
  }
}

/**
 * Returns the elementwise raising of the first argument to the power of the
 * second argument.
 *
 * @tparam T1 type of first argument
 * @tparam T2 type of second argument
 * @param a first argument
 * @param b second argument
 * @return the elementwise raising of the first argument to the power of the
 * second argument.
 */
template <typename T1, typename T2, require_any_container_t<T1, T2>* = nullptr,
          require_all_not_matrix_st<is_var, T1, T2>* = nullptr,
          require_any_fvar_t<base_type_t<T1>, base_type_t<T2>>* = nullptr>
inline auto pow(const T1& a, const T2& b) {
  return apply_scalar_binary(
      a, b, [](const auto& c, const auto& d) { return stan::math::pow(c, d); });
}

}  // namespace math
}  // namespace stan
#endif
