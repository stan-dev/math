#ifndef STAN_MATH_REV_FUN_BETA_HPP
#define STAN_MATH_REV_FUN_BETA_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/prim/fun/beta.hpp>
#include <stan/math/prim/fun/digamma.hpp>

namespace stan {
namespace math {

/**
 * Returns the beta function and gradients for two var inputs.
 *
   \f[
     \mathrm{beta}(a,b) = \left(B\left(a,b\right)\right)
   \f]

   \f[
    \frac{\partial }{\partial a} = \left(\psi^{\left(0\right)}\left(a\right)
                                      - \psi^{\left(0\right)}
                                      \left(a + b\right)\right)
                                    * \mathrm{beta}(a,b)
   \f]

   \f[
    \frac{\partial }{\partial b} = \left(\psi^{\left(0\right)}\left(b\right)
                                      - \psi^{\left(0\right)}
                                      \left(a + b\right)\right)
                                    * \mathrm{beta}(a,b)
   \f]
 *
 * @param a var Argument
 * @param b var Argument
 * @return Result of beta function
 */
inline var beta(const var& a, const var& b) {
  double digamma_ab = digamma(a.val() + b.val());
  double digamma_a = digamma(a.val()) - digamma_ab;
  double digamma_b = digamma(b.val()) - digamma_ab;
  return make_callback_var(beta(a.val(), b.val()),
                           [a, b, digamma_a, digamma_b](auto& vi) mutable {
                             const double adj_val = vi.adj() * vi.val();
                             a.adj() += adj_val * digamma_a;
                             b.adj() += adj_val * digamma_b;
                           });
}

/**
 * Returns the beta function and gradient for first var input.
 *
   \f[
     \mathrm{beta}(a,b) = \left(B\left(a,b\right)\right)
   \f]

   \f[
    \frac{\partial }{\partial a} = \left(\psi^{\left(0\right)}\left(a\right)
                                      - \psi^{\left(0\right)}
                                      \left(a + b\right)\right)
                                    * \mathrm{beta}(a,b)
   \f]
 *
 * @param a var Argument
 * @param b double Argument
 * @return Result of beta function
 */
inline var beta(const var& a, double b) {
  auto digamma_ab = digamma(a.val()) - digamma(a.val() + b);
  return make_callback_var(beta(a.val(), b),
                           [a, b, digamma_ab](auto& vi) mutable {
                             a.adj() += vi.adj() * digamma_ab * vi.val();
                           });
}

/**
 * Returns the beta function and gradient for second var input.
 *
   \f[
     \mathrm{beta}(a,b) = \left(B\left(a,b\right)\right)
   \f]

   \f[
    \frac{\partial }{\partial b} = \left(\psi^{\left(0\right)}\left(b\right)
                                      - \psi^{\left(0\right)}
                                      \left(a + b\right)\right)
                                    * \mathrm{beta}(a,b)
   \f]
 *
 * @param a double Argument
 * @param b var Argument
 * @return Result of beta function
 */
inline var beta(double a, const var& b) {
  auto beta_val = beta(a, b.val());
  auto digamma_ab = (digamma(b.val()) - digamma(a + b.val())) * beta_val;
  return make_callback_var(beta_val, [a, b, digamma_ab](auto& vi) mutable {
    b.adj() += vi.adj() * digamma_ab;
  });
}

}  // namespace math
}  // namespace stan
#endif
