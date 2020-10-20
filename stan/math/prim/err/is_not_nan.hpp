#ifndef STAN_MATH_PRIM_ERR_IS_NOT_NAN_HPP
#define STAN_MATH_PRIM_ERR_IS_NOT_NAN_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/get.hpp>
#include <stan/math/prim/fun/is_nan.hpp>
#include <stan/math/prim/fun/size.hpp>
#include <stan/math/prim/fun/value_of_rec.hpp>

namespace stan {
namespace math {

/**
 * Return <code>true</code> if <code>y</code> is not <code>NaN</code>.
 * This function is vectorized and will check each element of
 * <code>y</code>. If no element is <code>NaN</code>, this
 * function will return <code>true</code>.
 * @tparam T_y Type of y
 * @param y Variable to check
 * @return <code>true</code> if no element of y is NaN
 */
template <typename T_y>
inline bool is_not_nan(const T_y& y) {
  for (size_t n = 0; n < stan::math::size(y); ++n) {
    if (is_nan(value_of_rec(stan::get(y, n)))) {
      return false;
    }
  }
  return true;
}

}  // namespace math
}  // namespace stan
#endif
