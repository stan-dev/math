#ifndef STAN_MATH_PRIM_FUN_INV_INC_BETA_HPP
#define STAN_MATH_PRIM_FUN_INV_INC_BETA_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/boost_policy.hpp>
#include <boost/math/special_functions/beta.hpp>

namespace stan {
namespace math {

/**
 * The inverse of the normalized incomplete beta function of a, b, with
 * probability p.
 *
 * Used to compute the inverse cumulative density function for the beta
 * distribution.
 *
 * @param a Shape parameter a >= 0; a and b can't both be 0
 * @param b Shape parameter b >= 0
 * @param p Random variate. 0 <= p <= 1
 * @throws if constraints are violated or if any argument is NaN
 * @return The inverse of the normalized incomplete beta function.
 */
inline double inv_inc_beta(double a, double b, double p) {
  check_not_nan("inv_inc_beta", "a", a);
  check_not_nan("inv_inc_beta", "b", b);
  check_not_nan("inv_inc_beta", "p", p);
  check_positive("inv_inc_beta", "a", a);
  check_positive("inv_inc_beta", "b", b);
  check_bounded("inv_inc_beta", "p", p, 0, 1);
  return boost::math::ibeta_inv(a, b, p, boost_policy_t<>());
}

/**
 * Enables the vectorized application of the inv_inc_beta function, when
 *  any arguments are containers.
 *
 * @tparam T1 type of first input
 * @tparam T2 type of second input
 * @tparam T3 type of third input
 * @param a First input
 * @param b Second input
 * @param c Third input
 * @return Inverse of the incomplete Beta function applied to the three inputs.
 */
template <typename T1, typename T2, typename T3,
          require_any_container_t<T1, T2, T3>* = nullptr>
inline auto inv_inc_beta(const T1& a, const T2& b, const T3& c) {
  return apply_scalar_ternary(
      [](const auto& d, const auto& e, const auto& f) {
        return inv_inc_beta(d, e, f);
      },
      a, b, c);
}

}  // namespace math
}  // namespace stan
#endif
