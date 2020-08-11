#ifndef STAN_MATH_PRIM_ERR_CHECK_NOT_NAN_SCREEN_HPP
#define STAN_MATH_PRIM_ERR_CHECK_NOT_NAN_SCREEN_HPP

#include <stan/math/prim/err/elementwise_check.hpp>
#include <stan/math/prim/fun/value_of_rec.hpp>

namespace stan {
namespace math {

/**
 * Run a quick check on an eigen container to see if it contains
 *   any NaNs.
 *
 * If returns false, the values in the container are not NaNs.
 *
 * If it returns true, a more thorough check is necessary.
 *
 * @tparam T_y type of y
 * @param y variable to check
 * @return true if full check is necessary
 */
template <typename T_y, require_eigen_t<T_y>* = nullptr>
inline bool check_not_nan_screen(const T_y& y) {
  return value_of_rec(y).hasNaN();
}

/**
 * For input types where no quick NaN check is available, return true
 *  (implying the elements in the container could be NaN)
 * 
 * @tparam T_y type of y
 * @param y variable to check
 * @return true if full check is necessary
 */
template <typename T_y, require_not_eigen_t<T_y>* = nullptr>
inline bool check_not_nan_screen(const T_y& y) {
  return true;
}

}  // namespace math
}  // namespace stan

#endif
