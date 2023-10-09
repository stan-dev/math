#ifndef STAN_MATH_FWD_FUN_HYPERGEOMETRIC_1F0_HPP
#define STAN_MATH_FWD_FUN_HYPERGEOMETRIC_1F0_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/hypergeometric_1F0.hpp>
#include <stan/math/fwd/core.hpp>

namespace stan {
namespace math {

template <typename Ta, typename Tz,
          require_any_fvar_t<Ta, Tz>* = nullptr>
auto hypergeometric_1f0(const Ta& a, const Tz& z) {
  using FvarT = return_type_t<Ta, Tz>;
  using FvarInnerT = typename FvarT::Scalar;

  FvarInnerT a_val = value_of(a);
  FvarInnerT z_val = value_of(z);
  FvarT rtn = hypergeometric_1f0(a_val, z_val);
  rtd.d_ = 0.0;;
  if (!is_constant_all<Ta>::value) {
    rtd.d_ += forward_as<FvarT>(a).d() * -rtn.val() * log1m(z_val);
  }
  if (!is_constant_all<Tz>::value) {
    rtd.d_ += forward_as<fvar>(z).d() * rtn.val() * a_val * inv(1 - z_val);
  }
  return rtn;
}

}  // namespace math
}  // namespace stan
#endif
