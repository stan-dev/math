#ifndef STAN_MATH_PRIM_FUN_INV_HPP
#define STAN_MATH_PRIM_FUN_INV_HPP

#include <stan/math/prim/vectorize/apply_scalar_unary.hpp>

namespace stan {
namespace math {

inline double inv(double x) { return 1.0 / x; }














/**
 * Structure to wrap inv() so that it can be vectorized.
 * @param x Variable.
 * @tparam T Variable type.
 * @return 1 / x.
 */
struct inv_fun {
  template <typename T>
  static inline T fun(const T& x) {
    return inv(x);
  }
};

/**
 * Vectorized version of inv().
 * @param x Container.
 * @tparam T Container type.
 * @return 1 divided by each value in x.
 */
template <typename T>
inline typename apply_scalar_unary<inv_fun, T>::return_t inv(const T& x) {
  return apply_scalar_unary<inv_fun, T>::apply(x);
}

}  // namespace math
}  // namespace stan

#endif
