#ifndef STAN_MATH_PRIM_ERR_CHECK_FINITE_HPP
#define STAN_MATH_PRIM_ERR_CHECK_FINITE_HPP

#include <stan/math/prim/err/elementwise_check.hpp>
#include <stan/math/prim/err/check_finite_screen.hpp>

namespace stan {
namespace math {

/**
 * Check if <code>y</code> is finite.
 * This function is vectorized and will check each element of
 * <code>y</code>.
 * @tparam T_y type of y
 * @param function function name (for error messages)
 * @param name variable name (for error messages)
 * @param y variable to check
 * @throw <code>domain_error</code> if y is infinity, -infinity, or NaN
 */
template <typename T_y>
inline void check_finite(const char* function, const char* name, const T_y& y) {
  if (check_finite_screen(y)) {
    auto is_good = [](const auto& y) { return std::isfinite(y); };
    elementwise_check(is_good, function, name, y, ", but must be finite!");
  }
}

}  // namespace math
}  // namespace stan

#endif
