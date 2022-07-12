#ifndef STAN_MATH_PRIM_FUN_HYPERGEOMETRIC_2F2_HPP
#define STAN_MATH_PRIM_FUN_HYPERGEOMETRIC_2F2_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/hypergeometric_pFq.hpp>
#include <boost/math/special_functions/hypergeometric_pFq.hpp>

namespace stan {
namespace math {

/**
 * Returns the Gauss hypergeometric function applied to the
 * input arguments:
 * \f$_2F_1(a_1,a_2;b;z)\f$
 *
 * See 'grad_2F1.hpp' for the derivatives wrt each parameter
 *
 * @param[in] a1 First of 'a' arguments to function
 * @param[in] a2 Second of 'a' arguments to function
 * @param[in] b 'b' argument to function
 * @param[in] z Scalar z argument
 * @return Gauss hypergeometric function
 */
template <typename Ta1, typename Ta2, typename Tb, typename Tz,
          require_all_arithmetic_t<Ta1, Ta2, Tb, Tz>* = nullptr>
return_type_t<Ta1, Ta1, Tb, Tz> hypergeometric_2F1(const Ta1& a1, const Ta2& a2,
                                                   const Tb& b, const Tz& z) {
  check_finite("hypergeometric_2F1", "a1", a1);
  check_finite("hypergeometric_2F1", "a2", a2);
  check_finite("hypergeometric_2F1", "b", b);
  check_finite("hypergeometric_2F1", "z", z);

  check_not_nan("hypergeometric_2F1", "a1", a1);
  check_not_nan("hypergeometric_2F1", "a2", a2);
  check_not_nan("hypergeometric_2F1", "b", b);
  check_not_nan("hypergeometric_2F1", "z", z);

  check_2F1_converges("hypergeometric_2F1", a1, a2, b, z);

  return boost::math::hypergeometric_pFq({a1, a2}, {b}, z);
}
}  // namespace math
}  // namespace stan
#endif
