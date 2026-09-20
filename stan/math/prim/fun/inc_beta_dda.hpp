#ifndef STAN_MATH_PRIM_FUN_INC_BETA_DDA_HPP
#define STAN_MATH_PRIM_FUN_INC_BETA_DDA_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/digamma.hpp>
#include <stan/math/prim/fun/grad_reg_inc_beta.hpp>

namespace stan {
namespace math {

/**
 * Returns the partial derivative of the regularized
 * incomplete beta function, I_{z}(a, b) with respect to a.
 *
 * Evaluated by <code>grad_reg_inc_beta</code>, which computes both shape
 * derivatives; call that function directly when both are needed.
 *
 * @tparam T scalar types of arguments
 * @param a first argument
 * @param b second argument
 * @param z upper bound of the integral
 * @param digamma_a value of digamma(a)
 * @param digamma_ab value of digamma(a + b)
 * @return partial derivative of the incomplete beta with respect to a
 *
 * @pre a >= 0
 * @pre b >= 0
 * @pre 0 <= z <= 1
 */
template <typename T>
inline T inc_beta_dda(T a, T b, T z, T digamma_a, T digamma_ab) {
  T g1 = 0;
  T g2 = 0;
  grad_reg_inc_beta(g1, g2, a, b, z, digamma_a, digamma(b), digamma_ab, T(0));
  return g1;
}

}  // namespace math
}  // namespace stan
#endif
