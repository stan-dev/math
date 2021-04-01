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
 * @param a Variable argument.
 * @return The inverse complementary log-log of the specified
 * argument.
 */
inline var inv_cloglog(const var& a) {
  auto precomp_exp = std::exp(a.val() - std::exp(a.val()));
  return make_callback_var(inv_cloglog(a.val()),
                           [a, precomp_exp](auto& vi) mutable {
                             a.adj() += vi.adj() * precomp_exp;
                           });
}

}  // namespace math
}  // namespace stan
#endif
