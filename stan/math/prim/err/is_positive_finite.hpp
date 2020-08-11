#ifndef STAN_MATH_PRIM_ERR_IS_POSITIVE_FINITE_HPP
#define STAN_MATH_PRIM_ERR_IS_POSITIVE_FINITE_HPP

#include <stan/math/prim/err/elementwise_check.hpp>
#include <stan/math/prim/err/check_finite_screen.hpp>

namespace stan {
namespace math {

/**
 * Return <code>true</code> if <code>y</code> is positive and finite.
 * This function is vectorized and will check each element of
 * <code>y</code>.
 * @tparam T_y Type of y
 * @param y Variable to check
 * @return <code>true</code> if every element of y is >0 and if no element of y
 * is NaN or infinity.
 */
template <typename T_y>
inline bool is_positive_finite(const T_y& y) {
  if (check_finite_screen(y)) {
    auto is_good = [](const auto& y) { return y > 0 && std::isfinite(y); };
    return elementwise_is(is_good, y);
  }
}
}  // namespace math
}  // namespace stan
#endif
