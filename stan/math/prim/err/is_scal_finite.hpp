#ifndef STAN_MATH_PRIM_ERR_IS_SCAL_FINITE_HPP
#define STAN_MATH_PRIM_ERR_IS_SCAL_FINITE_HPP

#include <stan/math/prim/err/elementwise_check.hpp>
#include <stan/math/prim/err/check_finite_screen.hpp>

namespace stan {
namespace math {

/**
 * Return <code>true</code> if <code>y</code> is finite.
 * This function is vectorized and will check each element of
 * <code>y</code>.
 * @tparam T_y Type of y
 * @param y Variable to check
 * @return <code>true</code> if no element in y is infinity, -infinity, or NaN
 */
template <typename T_y>
inline bool is_scal_finite(const T_y& y) {
  if (check_finite_screen(y)) {
    auto is_good = [](const auto& y) { return std::isfinite(y); };
    return elementwise_is(is_good, y);
  }

  return true;
}

}  // namespace math
}  // namespace stan
#endif
