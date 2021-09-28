#ifndef STAN_MATH_PRIM_ERR_CHECK_NONNEGATIVE_HPP
#define STAN_MATH_PRIM_ERR_CHECK_NONNEGATIVE_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err/elementwise_check.hpp>
#include <stan/math/prim/fun/get.hpp>
#include <stan/math/prim/fun/size.hpp>
#include <type_traits>

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
  elementwise_check([](double x) { return x >= 0; }, function, name, y,
                    "nonnegative");
}
}  // namespace math
}  // namespace stan
#endif
