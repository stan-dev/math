#ifndef STAN_MATH_PRIM_FUN_HYPERGEOMETRIC_1F0_HPP
#define STAN_MATH_PRIM_FUN_HYPERGEOMETRIC_1F0_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err/check_less.hpp>
#include <stan/math/prim/fun/boost_policy.hpp>
#include <boost/math/special_functions/hypergeometric_1F0.hpp>
#include <cmath>

namespace stan {
namespace math {

/**
 * Returns the Hypergeometric 1F0 function applied to the
 * input arguments:
 * \f$ _1F_0(a;;z) = \sum_{k=1}^{\infty}\frac{\left(a\right)_kz^k}{k!}\f$
 *
 * \f$ \frac{\partial _1F_0\left(a;;z\right)}{\partial a} =
 *    -\left(1-z\right)^{-a}\log\left(1 - z\right) \f$
 *
 * \f$ \frac{\partial _1F_0\left(a;;z\right)}{\partial z} =
 *    a\left(1-z\right)^{-a-1} \f$
 *
 * @tparam Ta Arithmetic type of 'a' argument
 * @tparam Tz Arithmetic type of 'z' argument
 * @param[in] a Scalar 'a' argument
 * @param[in] z Scalar z argument
 * @return Hypergeometric 1F0 function
 */
template <typename Ta, typename Tz, require_all_arithmetic_t<Ta, Tz>* = nullptr>
return_type_t<Ta, Tz> hypergeometric_1f0(const Ta& a, const Tz& z) {
  constexpr const char* function = "hypergeometric_1f0";
  check_less("hypergeometric_1f0", "abs(z)", std::fabs(z), 1.0);

  return boost::math::hypergeometric_1F0(a, z, boost_policy_t<>());
}

}  // namespace math
}  // namespace stan
#endif
