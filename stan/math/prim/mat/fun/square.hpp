#ifndef STAN_MATH_PRIM_MAT_FUN_SQUARE_HPP
#define STAN_MATH_PRIM_MAT_FUN_SQUARE_HPP

#include <stan/math/prim/vectorize/apply_scalar_unary.hpp>
#include <stan/math/prim/scal/fun/square.hpp>

namespace stan {
namespace math {

/**
 * Structure to wrap square() so that it can be vectorized.
 *
 * @tparam T type of variable
 * @param x variable
 * @return x squared.
 */
struct square_fun {
  template <typename T>
  static inline T fun(const T& x) {
    return square(x);
  }
};

/**
 * Vectorized version of square().
 *
 * @tparam T type of container
 * @param x container
 * @return Each value in x squared.
 */
template <typename T>
inline typename apply_scalar_unary<square_fun, T>::return_t square(const T& x) {
  return apply_scalar_unary<square_fun, T>::apply(x);
}

}  // namespace math
}  // namespace stan

#endif
