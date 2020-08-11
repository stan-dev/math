#ifndef STAN_MATH_PRIM_ERR_CHECK_POSITIVE_FINITE_HPP
#define STAN_MATH_PRIM_ERR_CHECK_POSITIVE_FINITE_HPP

#include <stan/math/prim/err/elementwise_check.hpp>
#include <stan/math/prim/err/check_finite_screen.hpp>

namespace stan {
namespace math {

/**
 * Check if <code>y</code> is positive and finite.
 * This function is vectorized and will check each element of
 * <code>y</code>.
 * @tparam T_y Type of y
 * @param function Function name (for error messages)
 * @param name Variable name (for error messages)
 * @param y Variable to check
 * @throw <code>domain_error</code> if any element of y is not positive or
 *   if any element of y is NaN or infinity.
 */

template <typename T_y>
inline void check_positive_finite(const char* function, const char* name,
                                  const T_y& y) {
  if(check_finite_screen(y)) {
    auto is_good = [](const auto& y) { return y > 0 && std::isfinite(y); };
    elementwise_check(is_good, function, name, y,
		      ", but must be positive and finite!");
  }
}
}  // namespace math
}  // namespace stan
#endif
