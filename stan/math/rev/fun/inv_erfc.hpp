#ifndef STAN_MATH_REV_FUN_INV_ERFC_HPP
#define STAN_MATH_REV_FUN_INV_ERFC_HPP

#include <stan/math/rev/core.hpp>
#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/fun/square.hpp>
#include <stan/math/rev/fun/exp.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/inv_erfc.hpp>
#include <cmath>

namespace stan {
namespace math {

/**
 * The inverse complementary error function for variables.
 *
 * The derivative is:
 \f[
    \frac{d}{dx} \left(\mbox{erfc}^{-1}(x)\right) =
      -\frac{\sqrt{\pi}}{2}e^{\mbox{erfc}^{-1}(x)^2}
 \f].
 *
 *
 * @param a The variable.
 * @return Inverse complementary error function applied to the variable.
 */
inline var inv_erfc(const var& a) {
  auto precomp_inv_erfc = inv_erfc(a.val());
  return make_callback_var(
      precomp_inv_erfc, [a, precomp_inv_erfc](auto& vi) mutable {
        a.adj()
            -= vi.adj() * exp(LOG_SQRT_PI - LOG_TWO + square(precomp_inv_erfc));
      });
}

template <typename T, require_matrix_t<T>* = nullptr>
inline auto inv_erfc(const var_value<T>& a) {
  auto precomp_inv_erfc = to_arena(inv_erfc(a.val()));
  return make_callback_var(
      precomp_inv_erfc, [a, precomp_inv_erfc](auto& vi) mutable {
        a.adj().array()
            -= vi.adj().array()
               * exp(LOG_SQRT_PI - LOG_TWO + square(precomp_inv_erfc).array());
      });
}

}  // namespace math
}  // namespace stan
#endif
