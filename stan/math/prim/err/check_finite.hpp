#ifndef STAN_MATH_PRIM_ERR_CHECK_FINITE_HPP
#define STAN_MATH_PRIM_ERR_CHECK_FINITE_HPP

#include <stan/math/prim/err/elementwise_check.hpp>

namespace stan {
namespace math {
  
/**
 * Check if each element of the eigen variable <code>y</code> is finite.
 *
 * @tparam T_y type of y
 * @param function function name (for error messages)
 * @param name variable name (for error messages)
 * @param y variable to check
 * @throw <code>domain_error</code> if y is infinity, -infinity, or NaN
 */
template <typename T_y, require_eigen_t<T_y>* = nullptr>
inline void check_finite(const char* function, const char* name, const T_y& y) {
  if(!value_of_rec(y).allFinite()) {
    auto is_good = [](const auto& y) { return std::isfinite(y); };
    elementwise_check(is_good, function, name, y, ", but must be finite!");
  }
}

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
template <typename T_y, require_not_eigen_t<T_y>* = nullptr>
inline void check_finite(const char* function, const char* name, const T_y& y) {
  auto is_good = [](const auto& y) { return std::isfinite(y); };
  elementwise_check(is_good, function, name, y, ", but must be finite!");
}

}  // namespace math
}  // namespace stan

#endif
