#ifndef STAN_MATH_PRIM_FUN_HYPERGEOMETRIC_1F0_HPP
#define STAN_MATH_PRIM_FUN_HYPERGEOMETRIC_1F0_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/boost_policy.hpp>
#include <stan/math/prim/fun/floor.hpp>
#include <boost/math/special_functions/hypergeometric_1F0.hpp>

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
template <typename Ta, typename Tz,
          require_all_arithmetic_t<Ta, Tz>* = nullptr>
auto hypergeometric_1f0(const Ta& a, const Tz& z) {
  bool condition_1 = z == 1.0;
  bool condition_2 = (1.0 - z < 0.0) && floor(a) != a;
  if (condition_1 || condition_2) {
    std::stringstream msg;
    msg << "Hypergeometric 1F0 is undefined when z == 1.0 or "
        << "1 - z < 0 and a not an integer, but the following "
        << "arguments provided: "
        << "a: " << a << ", z: " << z
        << std::endl;
    throw std::domain_error(msg.str());
  }

  return boost::math::hypergeometric_1F0(a, z, boost_policy_t<>());
}

}  // namespace math
}  // namespace stan
#endif
