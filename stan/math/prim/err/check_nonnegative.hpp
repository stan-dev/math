#ifndef STAN_MATH_PRIM_ERR_CHECK_NONNEGATIVE_HPP
#define STAN_MATH_PRIM_ERR_CHECK_NONNEGATIVE_HPP

#include <stan/math/prim/err/elementwise_check.hpp>

namespace stan {
namespace math {

/**
 * Check if <code>y</code> is non-negative.
 * This function is vectorized and will check each element of <code>y</code>.
 * @tparam T_y Type of y
 * @param function Function name (for error messages)
 * @param name Variable name (for error messages)
 * @param y Variable to check
 * @throw <code>domain_error</code> if y is negative or
 *   if any element of y is NaN.
 */
template <typename T_y>
inline void check_nonnegative(const char* function, const char* name,
                              const T_y& y) {
  auto is_good = [](const auto& y) { return y >= 0; };
  elementwise_check(is_good, function, name, y, ", but must be >= 0!");
}

}  // namespace math
}  // namespace stan
#endif
