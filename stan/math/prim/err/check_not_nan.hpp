#ifndef STAN_MATH_PRIM_ERR_CHECK_NOT_NAN_HPP
#define STAN_MATH_PRIM_ERR_CHECK_NOT_NAN_HPP

#include <stan/math/prim/err/elementwise_check.hpp>
#include <stan/math/prim/err/check_not_nan_screen.hpp>

namespace stan {
namespace math {

/**
 * Check if <code>y</code> is not <code>NaN</code>.
 * This function is vectorized and will check each element of
 * <code>y</code>. If any element is <code>NaN</code>, this
 * function will throw an exception.
 * @tparam T_y Type of y
 * @param function Function name (for error messages)
 * @param name Variable name (for error messages)
 * @param y Variable to check
 * @throw <code>domain_error</code> if any element of y is NaN
 */
template <typename T_y>
inline void check_not_nan(const char* function, const char* name,
                          const T_y& y) {
  if (check_not_nan_screen(y)) {
    auto is_good = [](const auto& y) { return !std::isnan(y); };
    elementwise_check(is_good, function, name, y, ", but must not be nan!");
  }
}

}  // namespace math
}  // namespace stan
#endif
