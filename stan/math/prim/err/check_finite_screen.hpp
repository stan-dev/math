#ifndef STAN_MATH_PRIM_ERR_CHECK_FINITE_SCREEN_HPP
#define STAN_MATH_PRIM_ERR_CHECK_FINITE_SCREEN_HPP

#include <stan/math/prim/err/elementwise_check.hpp>
#include <stan/math/prim/fun/value_of_rec.hpp>

namespace stan {
namespace math {

/**
 * Run a quick check on an eigen container to see if the elements
 *  are finite.
 *
 * If returns false, the values in the container are finite.
 *
 * If it returns true, a more thorough check is necessary.
 *
 * @tparam T_y type of y
 * @param y variable to check
 * @return true if full check is necessary
 */
template <typename T_y, require_eigen_t<T_y>* = nullptr>
inline bool check_finite_screen(const T_y& y) {
  return !value_of_rec(y).allFinite();
}

/**
 * For input types where no quick finite check is available, return true
 *  (implying the elements in the container are not necessarily finite)
 *
 * @tparam T_y type of y
 * @param y variable to check
 * @return true if full check is necessary
 */
template <typename T_y, require_not_eigen_t<T_y>* = nullptr>
inline bool check_finite_screen(const T_y& y) {
  return true;
}

}  // namespace math
}  // namespace stan

#endif
