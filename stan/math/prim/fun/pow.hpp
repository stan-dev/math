#ifndef STAN_MATH_PRIM_FUN_POW_HPP
#define STAN_MATH_PRIM_FUN_POW_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/functor/apply_scalar_binary.hpp>
#include <cmath>
#include <complex>

namespace stan {
namespace math {

namespace internal {

/**
 * Return the first argument raised to the power of the second
 * argument.  At least one of the arguments must be a complex number.
 *
 * @tparam U type of base
 * @tparam V type of exponent
 * @param[in] x base
 * @param[in] y exponent
 * @return base raised to the power of the exponent
 */
template <typename U, typename V>
inline complex_return_t<U, V> complex_pow(const U& x, const V& y) {
  return exp(y * log(x));
}
}  // namespace internal

/**
 * Return the first argument raised to the power of the second
 * argument.
 *
 * @tparam T1 type of first argument
 * @tparam T2 type of second argument
 * @param a first argument
 * @param b second argument
 * @return the first argument raised to the power of the second
 * argument.
 */
template <typename T1, typename T2,
          require_all_t<
              disjunction<is_complex<T1>, std::is_arithmetic<T1>>,
              disjunction<is_complex<T2>, std::is_arithmetic<T2>>>* = nullptr>
inline auto pow(const T1& a, const T2& b) {
  return std::pow(a, b);
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
          require_all_not_matrix_st<is_var, T1, T2>* = nullptr>
inline auto pow(const T1& a, const T2& b) {
  return apply_scalar_binary(a, b, [](const auto& c, const auto& d) {
    using std::pow;
    return pow(c, d);
  });
}
}  // namespace math
}  // namespace stan
#endif
