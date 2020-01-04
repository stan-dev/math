#ifndef STAN_MATH_PRIM_MAT_FUN_FLOOR_HPP
#define STAN_MATH_PRIM_MAT_FUN_FLOOR_HPP

#include <stan/math/prim/vectorize/apply_scalar_unary.hpp>
#include <cmath>

namespace stan {
namespace math {

/**
 * Structure to wrap floor() so that it can be vectorized.
 *
 * @tparam T type of variable
 * @param x variable
 * @return Greatest integer <= x.
 */
struct floor_fun {
  template <typename T>
  static inline T fun(const T& x) {
    using std::floor;
    return floor(x);
  }
};

/**
 * Vectorized version of floor().
 *
 * @tparam T type of container
 * @param x container
 * @return Greatest integer <= each value in x.
 */
template <typename T>
inline typename apply_scalar_unary<floor_fun, T>::return_t floor(const T& x) {
  return apply_scalar_unary<floor_fun, T>::apply(x);
}

}  // namespace math
}  // namespace stan

#endif
