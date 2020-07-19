#ifndef STAN_MATH_PRIM_ERR_IS_POSITIVE_HPP
#define STAN_MATH_PRIM_ERR_IS_POSITIVE_HPP

#include <stan/math/prim/err/elementwise_check.hpp>

namespace stan {
namespace math {

/**
 * Return <code>true</code> if <code>y</code> is positive.
 * This function is vectorized and will check each element of
 * <code>y</code>.
 * @tparam T_y Type of y
 * @param y Variable to check
 * @return <code>true</code> if y contains only positive elements
 */
template <typename T_y>
inline bool is_positive(const T_y& y) {
  auto is_good = [](const auto& y) { return y > 0; };
  return elementwise_is(is_good, y);
}

/**
 * Return <code>true</code> if <code>size</code> is positive.
 * @param size Size value to check
 * @return <code>true</code> if <code>size</code> is not zero or negative
 */
inline bool is_positive(int size) { return size > 0; }

}  // namespace math
}  // namespace stan
#endif
