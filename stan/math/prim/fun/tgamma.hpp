#ifndef STAN_MATH_PRIM_FUN_TGAMMA_HPP
#define STAN_MATH_PRIM_FUN_TGAMMA_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/is_nonpositive_integer.hpp>
#include <stan/math/prim/functor/apply_scalar_unary.hpp>
#include <cmath>

namespace stan {
namespace math {

/**
 * Return the gamma function applied to the specified argument.
 *
 * @param x Argument.
 * @return The gamma function applied to argument.
 */
inline double tgamma(double x) {
  if (x == 0.0 || is_nonpositive_integer(x)) {
    throw_domain_error("tgamma", "x", x, "x == 0 or negative integer");
  }
  return std::tgamma(x);
}

/**
 * Structure to wrap tgamma() so that it can be vectorized.
 *
 * @tparam T type of variable
 * @param x variable
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
 *
 * @tparam T type of container
 * @param x container
 * @return Gamma function applied to each value in x.
 * @throw std::domain_error if any value is 0 or a negative integer
 */
template <
    typename T,
    require_all_not_nonscalar_prim_or_rev_kernel_expression_t<T>* = nullptr,
    require_not_var_matrix_t<T>* = nullptr>
inline auto tgamma(const T& x) {
  return apply_scalar_unary<tgamma_fun, T>::apply(x);
}

}  // namespace math
}  // namespace stan

#endif
