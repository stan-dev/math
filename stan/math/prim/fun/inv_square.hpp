#ifndef STAN_MATH_PRIM_FUN_INV_SQUARE_HPP
#define STAN_MATH_PRIM_FUN_INV_SQUARE_HPP

#include <stan/math/prim/vectorize/apply_scalar_unary.hpp>

namespace stan {
namespace math {

inline double inv_square(double x) { return 1.0 / (x * x); }

/**
 * Structure to wrap inv_square() so that it can be vectorized.
 * @param x Variable.
 * @tparam T Variable type.
 * @return 1 / x squared.
 */
struct inv_square_fun {
  template <typename T>
  static inline T fun(const T& x) {
    return inv_square(x);
  }
};

/**
 * Vectorized version of inv_square().
 * @param x Container.
 * @tparam T Container type.
 * @return 1 / the square of each value in x.
 */
template <typename T>
inline typename apply_scalar_unary<inv_square_fun, T>::return_t inv_square(
    const T& x) {
  return apply_scalar_unary<inv_square_fun, T>::apply(x);
}

}  // namespace math
}  // namespace stan

#endif
