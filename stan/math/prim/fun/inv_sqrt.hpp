#ifndef STAN_MATH_PRIM_FUN_INV_SQRT_HPP
#define STAN_MATH_PRIM_FUN_INV_SQRT_HPP

#include <cmath>
#include <stan/math/prim/vectorize/apply_scalar_unary.hpp>

namespace stan {
namespace math {

inline double inv_sqrt(double x) {
  using std::sqrt;
  return 1.0 / sqrt(x);
}

/**
 * Structure to wrap inv_sqrt() so that it can be vectorized.
 * @param x Variable.
 * @tparam T Variable type.
 * @return 1 / sqrt of x.
 */
struct inv_sqrt_fun {
  template <typename T>
  static inline T fun(const T& x) {
    return inv_sqrt(x);
  }
};

/**
 * Vectorized version of inv_sqrt().
 * @param x Container.
 * @tparam T Container type.
 * @return 1 / sqrt of each value in x.
 */
template <typename T>
inline typename apply_scalar_unary<inv_sqrt_fun, T>::return_t inv_sqrt(
    const T& x) {
  return apply_scalar_unary<inv_sqrt_fun, T>::apply(x);
}

}  // namespace math
}  // namespace stan

#endif
