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
  return make_callback_var(beta(a.val(), b.val()), [a, b](auto& vi) mutable {
    const double adj_val = vi.adj() * vi.val();
    const double digamma_ab = digamma(a.val() + b.val());
    a.adj() += adj_val * (digamma(a.val()) - digamma_ab);
    b.adj() += adj_val * (digamma(b.val()) - digamma_ab);
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
  return make_callback_var(beta(a.val(), b), [a, b](auto& vi) mutable {
    a.adj() += vi.adj() * (digamma(a.val()) - digamma(a.val() + b)) * vi.val();
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
  return make_callback_var(beta(a, b.val()), [a, b](auto& vi) mutable {
    b.adj() += vi.adj() * (digamma(b.val()) - digamma(a + b.val())) * vi.val();
  });
}

}  // namespace math
}  // namespace stan
#endif
