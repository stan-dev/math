#ifndef STAN_MATH_PRIM_ERR_IS_NONNEGATIVE_HPP
#define STAN_MATH_PRIM_ERR_IS_NONNEGATIVE_HPP

#include <stan/math/prim/err/elementwise_check.hpp>

namespace stan {
namespace math {

/**
 * Return <code>true</code> if <code>y</code> is nonnegative.
 * This function is vectorized and will check each element of
 * <code>y</code>.
 * @tparam T_y Type of y
 * @param y Variable to check
 * @return <code>true</code> if every element of y is >=0.
 */
template <typename T_y>
inline bool is_nonnegative(const T_y& y) {
  auto is_good = [](const auto& y) { return y >= 0; };
  return elementwise_is(is_good, y);
}
}  // namespace math
}  // namespace stan
#endif
