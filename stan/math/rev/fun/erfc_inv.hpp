#ifndef STAN_MATH_REV_FUN_ERFC_INV_HPP
#define STAN_MATH_REV_FUN_ERFC_INV_HPP

#include <stan/math/rev/core.hpp>
#include <stan/math/rev/meta.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/erfc_inv.hpp>
#include <stan/math/prim/fun/square.hpp>
#include <cmath>

namespace stan {
namespace math {

/**
 * The inverse complementary error function for variables (C99).
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
inline var erfc_inv(const var& a) {
  auto precomp_erfc_inv = erfc_inv(a.val());
  return make_callback_var(precomp_erfc_inv,
                           [a, precomp_erfc_inv](auto& vi) mutable {
    a.adj() -= vi.adj()
                * exp(LOG_SQRT_PI - LOG_TWO + square(precomp_erfc_inv));
  });
}

template <typename T, require_matrix_t<T>* = nullptr>
inline auto erfc_inv(const var_value<T>& a) {
  auto precomp_erfc_inv = to_arena(erfc_inv(a.val()));
  return make_callback_var(precomp_erfc_inv,
                           [a, precomp_erfc_inv](auto& vi) mutable {
    a.adj().array() -= vi.adj().array()
            * exp(LOG_SQRT_PI - LOG_TWO + square(precomp_erfc_inv).array());
  });
}

}  // namespace math
}  // namespace stan
#endif
