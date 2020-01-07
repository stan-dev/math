#ifndef STAN_MATH_PRIM_MAT_FUN_CBRT_HPP
#define STAN_MATH_PRIM_MAT_FUN_CBRT_HPP

#include <stan/math/prim/vectorize/apply_scalar_unary.hpp>
#include <stan/math/prim/scal/fun/cbrt.hpp>

namespace stan {
namespace math {

/**
 * Structure to wrap cbrt() so it can be vectorized.
 *
 * @tparam T type of variable
 * @param x variable
 * @return Cube root of x.
 */
struct cbrt_fun {
  template <typename T>
  static inline T fun(const T& x) {
    return cbrt(x);
  }
};

/**
 * Vectorized version of cbrt().
 *
 * @tparam T type of container
 * @param x container
 * @return Cube root of each value in x.
 */
template <typename T>
inline typename apply_scalar_unary<cbrt_fun, T>::return_t cbrt(const T& x) {
  return apply_scalar_unary<cbrt_fun, T>::apply(x);
}

}  // namespace math
}  // namespace stan

#endif
