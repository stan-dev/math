#ifndef STAN_MATH_PRIM_SCAL_FUN_LDEXP_HPP
#define STAN_MATH_PRIM_SCAL_FUN_LDEXP_HPP

#include <stan/math/prim/meta.hpp>
#include <cmath>

namespace stan {
namespace math {

/**
 * Returns the product of a (the significand) and
 * 2 to power b (the exponent).
 *
 * @tparam T Scalar type of significand
 * @param[in] a the significand
 * @param[in] b an integer that is the exponent
 * @return product of a times 2 to the power b
 */
template <typename T1, typename T2,
          typename = enable_if_all_arithmetic<T1, T2>>
inline auto ldexp(T1&& a, T2&& b) {
  using std::ldexp;
  return ldexp(std::forward<T1>(a), std::forward<T2>(b));
}
}  // namespace math
}  // namespace stan

#endif
