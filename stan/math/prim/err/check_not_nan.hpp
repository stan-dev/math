#ifndef STAN_MATH_PRIM_ERR_CHECK_NOT_NAN_HPP
#define STAN_MATH_PRIM_ERR_CHECK_NOT_NAN_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err/elementwise_check.hpp>
#include <stan/math/prim/fun/get.hpp>
#include <stan/math/prim/fun/is_nan.hpp>
#include <stan/math/prim/fun/size.hpp>
#include <stan/math/prim/fun/value_of_rec.hpp>

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
  elementwise_check([](double x) { return !std::isnan(x); }, function, name, y,
                    "not nan");
}

}  // namespace math
}  // namespace stan
#endif
