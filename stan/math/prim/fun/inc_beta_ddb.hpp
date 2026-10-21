#ifndef STAN_MATH_PRIM_FUN_INC_BETA_DDB_HPP
#define STAN_MATH_PRIM_FUN_INC_BETA_DDB_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/grad_reg_inc_beta.hpp>

namespace stan {
namespace math {

/**
 * Returns the partial derivative of the regularized
 * incomplete beta function, I_{z}(a, b) with respect to b.
 *
 * Evaluated by <code>grad_reg_inc_beta</code>, computing only the derivative
 * with respect to b.
 *
 * @tparam T1 type of a
 * @tparam T2 type of b
 * @tparam T3 type of z
 * @tparam T4 type of digamma_b
 * @tparam T5 type of digamma_ab
 * @param a first argument
 * @param b second argument
 * @param z upper bound of the integral
 * @param digamma_b value of digamma(b)
 * @param digamma_ab value of digamma(a + b)
 * @return partial derivative of the incomplete beta with respect to b
 *
 * @pre a >= 0
 * @pre b >= 0
 * @pre 0 <= z <= 1
 */
template <typename T1, typename T2, typename T3, typename T4, typename T5,
          require_all_stan_scalar_t<T1, T2, T3, T4, T5>* = nullptr>
inline auto inc_beta_ddb(const T1& a, const T2& b, const T3& z,
                         const T4& digamma_b, const T5& digamma_ab) {
  return grad_reg_inc_beta<GRAD_RIB::G2>(a, b, z, 0.0, digamma_b, digamma_ab,
                                         0.0);
}

}  // namespace math
}  // namespace stan
#endif
