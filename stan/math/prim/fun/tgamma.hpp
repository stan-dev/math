#ifndef STAN_MATH_PRIM_FUN_TGAMMA_HPP
#define STAN_MATH_PRIM_FUN_TGAMMA_HPP

#include <stan/math/prim/err/domain_error.hpp>
#include <stan/math/prim/fun/is_nonpositive_integer.hpp>
#include <cmath>
#include <stan/math/prim/vectorize/apply_scalar_unary.hpp>
#include <stan/math/prim/fun/tgamma.hpp>

namespace stan {
namespace math {

/**
 * Return the gamma function applied to the specified argument.
 *
 * @param x Argument.
 * @return The gamma function applied to argument.
 */
inline double tgamma(double x) {
  if (x == 0.0 || is_nonpositive_integer(x))
    domain_error("tgamma", "x", x, "x == 0 or negative integer");
  return std::tgamma(x);
}

}  // namespace math
}  // namespace stan

namespace stan {
namespace math {

/**
 * Structure to wrap tgamma() so that it can be vectorized.
 * @param x Variable.
 * @tparam T Variable type.
 * @return Gamma function applied to x.
 * @throw std::domain_error if x is 0 or a negative integer
 */
struct tgamma_fun {
  template <typename T>
  static inline T fun(const T& x) {
    return tgamma(x);
  }
};

/**
 * Vectorized version of tgamma().
 * @param x Container.
 * @tparam T Container type.
 * @return Gamma function applied to each value in x.
 * @throw std::domain_error if any value is 0 or a negative integer
 */
template <typename T>
inline typename apply_scalar_unary<tgamma_fun, T>::return_t tgamma(const T& x) {
  return apply_scalar_unary<tgamma_fun, T>::apply(x);
}

}  // namespace math
}  // namespace stan

#endif
