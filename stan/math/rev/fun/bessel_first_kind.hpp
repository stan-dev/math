#ifndef STAN_MATH_REV_FUN_BESSEL_FIRST_KIND_HPP
#define STAN_MATH_REV_FUN_BESSEL_FIRST_KIND_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/prim/fun/bessel_first_kind.hpp>

namespace stan {
namespace math {

inline var bessel_first_kind(int v, const var& a) {
  double ret_val = bessel_first_kind(v, a.val());
  auto precomp_bessel = v * ret_val / a.val()
                        - bessel_first_kind(v + 1, a.val());
  return make_callback_var(ret_val,
                           [precomp_bessel, a](const auto& vi) mutable {
                             a.adj() += vi.adj_ * precomp_bessel;
                           });
}

}  // namespace math
}  // namespace stan
#endif
