#ifndef STAN_MATH_FWD_FUN_HYPERGEOMETRIC_1F0_HPP
#define STAN_MATH_FWD_FUN_HYPERGEOMETRIC_1F0_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/hypergeometric_1F0.hpp>
#include <stan/math/fwd/core.hpp>

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
 * @tparam Ta Fvar or arithmetic type of 'a' argument
 * @tparam Tz Fvar or arithmetic type of 'z' argument
 * @param[in] a Scalar 'a' argument
 * @param[in] z Scalar z argument
 * @return Hypergeometric 1F0 function
 */
template <typename Ta, typename Tz, typename FvarT = return_type_t<Ta, Tz>,
          require_all_stan_scalar_t<Ta, Tz>* = nullptr,
          require_any_fvar_t<Ta, Tz>* = nullptr>
FvarT hypergeometric_1f0(const Ta& a, const Tz& z) {
  partials_type_t<Ta> a_val = value_of(a);
  partials_type_t<Tz> z_val = value_of(z);
  FvarT rtn = FvarT(hypergeometric_1f0(a_val, z_val), 0.0);
  if (!is_constant_all<Ta>::value) {
    rtn.d_ += forward_as<FvarT>(a).d() * -rtn.val() * log1m(z_val);
  }
  if (!is_constant_all<Tz>::value) {
    rtn.d_ += forward_as<FvarT>(z).d() * rtn.val() * a_val * inv(1 - z_val);
  }
  return rtn;
}

}  // namespace math
}  // namespace stan
#endif
