#ifndef STAN_MATH_PRIM_FUN_LDEXP_HPP
#define STAN_MATH_PRIM_FUN_LDEXP_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/functor/apply_scalar_binary.hpp>
#include <cmath>

namespace stan {
namespace math {

/**
 * Returns the product of a (the significand) and
 * 2 to power b (the exponent).
 *
 * @tparam T1 scalar type of significand
 * @param[in] a the significand
 * @param[in] b an integer that is the exponent
 * @return product of a times 2 to the power b
 */
template <typename T1, require_arithmetic_t<T1>* = nullptr>
inline double ldexp(T1 a, int b) {
  using std::ldexp;
  return ldexp(a, b);
}

/**
 * Enables the vectorised application of the ldexp function,
 * when the first and/or second arguments are containers.
 *
 * @tparam T1 type of first input
 * @tparam T2 type of second input
 * @param a First input
 * @param b Second input
 * @return ldexp function applied to the two inputs.
 */
template <typename T1, typename T2, require_any_container_t<T1, T2>* = nullptr,
          require_all_not_nonscalar_prim_or_rev_kernel_expression_t<
              T1, T2>* = nullptr>
inline auto ldexp(const T1& a, const T2& b) {
  return apply_scalar_binary(
      a, b, [&](const auto& c, const auto& d) { return ldexp(c, d); });
}

}  // namespace math
}  // namespace stan

#endif
