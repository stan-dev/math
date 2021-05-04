#ifndef STAN_MATH_PRIM_ERR_CHECK_POSITIVE_FINITE_HPP
#define STAN_MATH_PRIM_ERR_CHECK_POSITIVE_FINITE_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err/elementwise_check.hpp>

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
 *   if any element of y is NaN.
 */
template <typename T_y>
inline void check_positive_finite(const char* function, const char* name,
                                  const T_y& y) {
  elementwise_check([](double x) { return x > 0 && std::isfinite(x); },
                    function, name, y, "positive finite");
}

}  // namespace math
}  // namespace stan
#endif
