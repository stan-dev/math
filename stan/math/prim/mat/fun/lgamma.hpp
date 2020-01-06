#ifndef STAN_MATH_PRIM_MAT_FUN_LGAMMA_HPP
#define STAN_MATH_PRIM_MAT_FUN_LGAMMA_HPP

#include <stan/math/prim/vectorize/apply_scalar_unary.hpp>
#include <stan/math/prim/scal/fun/lgamma.hpp>

namespace stan {
namespace math {

/**
 * Structure to wrap lgamma() so that it can be vectorized.
 *
 * @tparam T type of variable
 * @param x variable
 * @return Natural log of the gamma function applied to x.
 * @throw std::domain_error if x is a negative integer or 0.
 */
struct lgamma_fun {
  template <typename T>
  static inline T fun(const T& x) {
    return lgamma(x);
  }
};

/**
 * Vectorized version of lgamma().
 *
 * @tparam T type of container
 * @param x container
 * @return Natural log of the gamma function
 *         applied to each value in x.
 * @throw std::domain_error if any value is a negative integer or 0.
 */
template <typename T>
inline typename apply_scalar_unary<lgamma_fun, T>::return_t lgamma(const T& x) {
  return apply_scalar_unary<lgamma_fun, T>::apply(x);
}

}  // namespace math
}  // namespace stan

#endif
