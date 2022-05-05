#ifndef STAN_MATH_REV_FUN_INV_CLOGLOG_HPP
#define STAN_MATH_REV_FUN_INV_CLOGLOG_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/prim/fun/inv_cloglog.hpp>
#include <cmath>

namespace stan {
namespace math {

/**
 * Return the inverse complementary log-log function applied
 * specified variable (stan).
 *
 * See inv_cloglog() for the double-based version.
 *
 * The derivative is given by
 *
 * \f$\frac{d}{dx} \mbox{cloglog}^{-1}(x) = \exp (x - \exp (x))\f$.
 *
 * @tparam T Arithmetic or a type inheriting from `EigenBase`.
 * @param a Variable argument.
 * @return The inverse complementary log-log of the specified
 * argument.
 */
inline auto inv_cloglog(const var a) {
  auto precomp_exp = exp(a.val() - exp(a.val()));
  return make_callback_var(inv_cloglog(a.val()),
                           [a, precomp_exp](auto& vi) mutable {
                             a.adj() += vi.adj() * precomp_exp;
                           });
}

template <typename T, require_rev_matrix_t<T>* = nullptr>
inline auto inv_cloglog(const T& x) {
  auto x_arena = to_arena(x);
  auto precomp_exp
      = to_arena((exp(x_arena.val() - exp(x_arena.val()))).array());
  return make_callback_rev_matrix<T>(
      inv_cloglog(x_arena.val()), [x_arena, precomp_exp](auto& vi) mutable {
        x_arena.adj().array() += vi.adj().array() * precomp_exp;
      });
}

}  // namespace math
}  // namespace stan
#endif
