#ifndef STAN_MATH_PRIM_FUN_HYPERGEOMETRIC_2F2_HPP
#define STAN_MATH_PRIM_FUN_HYPERGEOMETRIC_2F2_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/hypergeometric_pFq.hpp>

namespace stan {
namespace math {

/**
 * Returns the generalised hypergeometric function applied to the
 * input arguments:
 * \f$_2F_2(a_1,a_2;b_1,b_2;z)\f$
 *
 * See 'grad_pFq.hpp' for the derivatives wrt each parameter
 *
 * @param[in] a Vector of 'a' arguments to function
 * @param[in] b Vector of 'b' arguments to function
 * @param[in] z Scalar z argument
 * @return Generalised hypergeometric function
 */
template <typename Ta, typename Tb, typename Tz,
          require_all_eigen_t<Ta, Tb>* = nullptr,
          require_stan_scalar_t<Tz>* = nullptr>
double hypergeometric_2F2(const Ta& a, const Tb& b, const Tz& z) {
   return hypergeometric_pFq(a, b, z);
}
}  // namespace math
}  // namespace stan
#endif
