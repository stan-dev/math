#ifndef STAN_MATH_PRIM_MAT_FUN_LAMBERTW_HPP
#define STAN_MATH_PRIM_MAT_FUN_LAMBERTW_HPP

#include <stan/math/prim/mat/vectorize/apply_scalar_unary.hpp>
#include <stan/math/scal.hpp>
#include <cmath>

namespace stan {
namespace math {

/**
 * Structure to wrap lambert_w0() so it can be vectorized.
 * @param x Variable.
 * @tparam T Variable type.
 * @return
 */
struct lam0_fun {
  template <typename T>
  static inline T fun(const T& x) {
    return lambert_w0(x);
  }
};

struct lam1_fun {
  template <typename T>
  static inline T fun(const T& x) {
    return lambert_wm1(x);
  }
};

/**
 * Vectorized version of lambert_w0().
 * @param x Container.
 * @tparam T Container type.
 * @return
 */
template <typename T>
inline typename apply_scalar_unary<lam0_fun, T>::return_t lambert_w0(const T& x) {
  return apply_scalar_unary<lam0_fun, T>::apply(x);
}

/**
 * Vectorized version of lambert_wm1().
 * @param x Container.
 * @tparam T Container type.
 * @return
 */
template <typename T>
inline typename apply_scalar_unary<lam1_fun, T>::return_t lambert_wm1(const T& x) {
  return apply_scalar_unary<lam1_fun, T>::apply(x);
}

}  // namespace math
}  // namespace stan

#endif
