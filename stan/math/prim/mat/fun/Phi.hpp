#ifndef STAN_MATH_PRIM_MAT_FUN_PHI_HPP
#define STAN_MATH_PRIM_MAT_FUN_PHI_HPP

#include <stan/math/prim/vectorize/apply_scalar_unary.hpp>
#include <stan/math/prim/scal/fun/Phi.hpp>

namespace stan {
namespace math {

/**
 * Structure to wrap Phi() so it can be vectorized.
 *
 * @tparam T type of argument
 * @param x argument
 * @return Unit normal CDF of x.
 */
struct Phi_fun {
  template <typename T>
  static inline T fun(const T& x) {
    return Phi(x);
  }
};

/**
 * Vectorized version of Phi().
 *
 * @tparam T type of container
 * @param x container
 * @return Unit normal CDF of each value in x.
 */
template <typename T>
inline typename apply_scalar_unary<Phi_fun, T>::return_t Phi(const T& x) {
  return apply_scalar_unary<Phi_fun, T>::apply(x);
}

}  // namespace math
}  // namespace stan

#endif
