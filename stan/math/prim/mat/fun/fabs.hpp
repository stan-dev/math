#ifndef STAN_MATH_PRIM_MAT_FUN_FABS_HPP
#define STAN_MATH_PRIM_MAT_FUN_FABS_HPP

#include <stan/math/prim/vectorize/apply_scalar_unary.hpp>
#include <cmath>

namespace stan {
namespace math {

/**
 * Structure to wrap fabs() so that it can be vectorized.
 *
 * @tparam T type of variable
 * @param x variable
 * @return Absolute value of x.
 */
struct fabs_fun {
  template <typename T>
  static inline T fun(const T& x) {
    using std::fabs;
    return fabs(x);
  }
};

/**
 * Vectorized version of fabs().
 *
 * @tparam T type of container
 * @param x container
 * @return Absolute value of each value in x.
 */
template <typename T>
inline typename apply_scalar_unary<fabs_fun, T>::return_t fabs(const T& x) {
  return apply_scalar_unary<fabs_fun, T>::apply(x);
}

}  // namespace math
}  // namespace stan

#endif
