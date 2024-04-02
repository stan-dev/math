#ifndef STAN_MATH_REV_FUN_HYPERGEOMETRIC_1F0_HPP
#define STAN_MATH_REV_FUN_HYPERGEOMETRIC_1F0_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/hypergeometric_1F0.hpp>
#include <stan/math/prim/fun/value_of.hpp>
#include <stan/math/prim/fun/log1m.hpp>
#include <stan/math/prim/fun/inv.hpp>
#include <stan/math/rev/core.hpp>

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
 * @tparam Ta Var or arithmetic type of 'a' argument
 * @tparam Tz Var or arithmetic type of 'z' argument
 * @param[in] a Scalar 'a' argument
 * @param[in] z Scalar z argument
 * @return Hypergeometric 1F0 function
 */
template <typename Ta, typename Tz,
          require_all_stan_scalar_t<Ta, Tz>* = nullptr,
          require_any_var_t<Ta, Tz>* = nullptr>
var hypergeometric_1f0(const Ta& a, const Tz& z) {
  double a_val = value_of(a);
  double z_val = value_of(z);
  double rtn = hypergeometric_1f0(a_val, z_val);
  return make_callback_var(rtn, [rtn, a, z, a_val, z_val](auto& vi) mutable {
    if (!is_constant_all<Ta>::value) {
      forward_as<var>(a).adj() += vi.adj() * -rtn * log1m(z_val);
    }
    if (!is_constant_all<Tz>::value) {
      forward_as<var>(z).adj() += vi.adj() * rtn * a_val * inv(1 - z_val);
    }
  });
}

}  // namespace math
}  // namespace stan
#endif
