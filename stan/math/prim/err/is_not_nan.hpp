#ifndef STAN_MATH_PRIM_ERR_IS_NOT_NAN_HPP
#define STAN_MATH_PRIM_ERR_IS_NOT_NAN_HPP

#include <stan/math/prim/err/elementwise_check.hpp>
#include <stan/math/prim/err/check_not_nan_screen.hpp>

namespace stan {
namespace math {

/**
 * Return <code>true</code> if <code>y</code> is not <code>NaN</code>.
 * This function is vectorized and will check each element of
 * <code>y</code>. If no element is <code>NaN</code>, this
 * function will return <code>true</code>.
 * @tparam T_y type of y
 * @param y variable to check
 * @return <code>true</code> if no element of y is NaN
 */
template <typename T_y>
inline bool is_not_nan(const T_y& y) {
  if (check_not_nan_screen(y)) {
    auto is_good = [](const auto& y) { return !std::isnan(y); };
    return elementwise_is(is_good, y);
  }
}

}  // namespace math
}  // namespace stan
#endif
