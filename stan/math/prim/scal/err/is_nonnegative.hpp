#ifndef STAN_MATH_PRIM_SCAL_ERR_IS_NONNEGATIVE_HPP
#define STAN_MATH_PRIM_SCAL_ERR_IS_NONNEGATIVE_HPP

#include <stan/math/prim/scal/meta/get.hpp>
#include <stan/math/prim/scal/meta/length.hpp>
#include <stan/math/prim/scal/meta/value_type.hpp>
#include <stan/math/prim/scal/meta/is_vector_like.hpp>
#include <type_traits>

namespace stan {
namespace math {

/**
 * Check if <code>y</code> is non-negative. 
 * This function is vectorized and will check each element 
 * of <code>y</code>.
 * @tparam T_y Type of <code>y</code>
 * @param y Variable to check
 * @return <code>true</code> if <code>y</code> is not negative or if no element
 *   of <code>y</code> is NaN
 */
template <typename T_y>
inline bool is_nonnegative(const T_y& y) {
  if (!(y >= 0))
    return false;

  for (size_t n = 0; n < stan::length(y); n++) {
    if (!(stan::get(y, n) >= 0))
      return false;
  }
  return true;
}

}  // namespace math
}  // namespace stan
#endif
