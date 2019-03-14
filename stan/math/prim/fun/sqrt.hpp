#ifndef STAN_MATH_PRIM_FUN_SQRT_HPP
#define STAN_MATH_PRIM_FUN_SQRT_HPP

#include <cmath>
#include <stan/math/prim/vectorize/apply_scalar_unary.hpp>



namespace stan {
namespace math {

/**
 * Return the square root of the specified argument.  This
 * version is required to disambiguate <code>sqrt(int)</code>.
 *
 * @param[in] x Argument.
 * @return Natural exponential of argument.
 */
inline double sqrt(int x) { return std::sqrt(x); }













/**
 * Structure to wrap sqrt() so that it can be vectorized.
 * @param x Variable.
 * @tparam T Variable type.
 * @return Square root of x.
 */
struct sqrt_fun {
  template <typename T>
  static inline T fun(const T& x) {
    using std::sqrt;
    return sqrt(x);
  }
};

/**
 * Vectorized version of sqrt().
 * @param x Container.
 * @tparam T Container type.
 * @return Square root of each value in x.
 */
template <typename T>
inline typename apply_scalar_unary<sqrt_fun, T>::return_t sqrt(const T& x) {
  return apply_scalar_unary<sqrt_fun, T>::apply(x);
}

}  // namespace math
}  // namespace stan

#endif
