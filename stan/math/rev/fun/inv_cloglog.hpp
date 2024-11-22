#ifndef STAN_MATH_REV_FUN_INV_CLOGLOG_HPP
#define STAN_MATH_REV_FUN_INV_CLOGLOG_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/fun/exp.hpp>
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
template <typename T, require_stan_scalar_or_eigen_t<T>* = nullptr>
inline auto inv_cloglog(const var_value<T>& a) {
  auto precomp_exp = to_arena(as_array_or_scalar(exp(a.val() - exp(a.val()))));
  return make_callback_var(inv_cloglog(a.val()),
                           [a, precomp_exp](auto& vi) mutable {
                             as_array_or_scalar(a.adj())
                                 += as_array_or_scalar(vi.adj()) * precomp_exp;
                           });
}

}  // namespace math
}  // namespace stan
#endif
