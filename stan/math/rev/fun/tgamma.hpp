#ifndef STAN_MATH_REV_FUN_TGAMMA_HPP
#define STAN_MATH_REV_FUN_TGAMMA_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/fun/digamma.hpp>
#include <stan/math/prim/fun/tgamma.hpp>

namespace stan {
namespace math {

/**
 * Return the Gamma function applied to the specified variable (C99).
 *
 * The derivative with respect to the argument is
 *
 * \f$\frac{d}{dx} \Gamma(x) = \Gamma(x) \Psi^{(0)}(x)\f$
 *
 * where \f$\Psi^{(0)}(x)\f$ is the digamma function.
 *
   \f[
   \mbox{tgamma}(x) =
   \begin{cases}
     \textrm{error} & \mbox{if } x\in \{\dots, -3, -2, -1, 0\}\\
     \Gamma(x) & \mbox{if } x\not\in \{\dots, -3, -2, -1, 0\}\\[6pt]
     \textrm{NaN} & \mbox{if } x = \textrm{NaN}
   \end{cases}
   \f]

   \f[
   \frac{\partial\, \mbox{tgamma}(x)}{\partial x} =
   \begin{cases}
     \textrm{error} & \mbox{if } x\in \{\dots, -3, -2, -1, 0\}\\
     \frac{\partial\, \Gamma(x)}{\partial x} & \mbox{if } x\not\in \{\dots, -3,
 -2, -1, 0\}\\[6pt] \textrm{NaN} & \mbox{if } x = \textrm{NaN} \end{cases} \f]

   \f[
   \Gamma(x)=\int_0^{\infty} u^{x - 1} \exp(-u) \, du
   \f]

   \f[
   \frac{\partial \, \Gamma(x)}{\partial x} = \Gamma(x)\Psi(x)
   \f]
 *
 * @param a Argument to function.
 * @return The Gamma function applied to the specified argument.
 */
inline var tgamma(const var& a) {
  return make_callback_var(tgamma(a.val()), [a](auto& vi) mutable {
    a.adj() += vi.adj() * vi.val() * digamma(a.val());
  });
}

/**
 * Return elementwise gamma function
 *
 * @tparam T a `var_value` with inner Eigen type
 * @param a input
 * @return elementwise gamma
 */
template <typename T, require_var_matrix_t<T>* = nullptr>
inline auto tgamma(const T& a) {
  return make_callback_var(tgamma(a.val()), [a](auto& vi) mutable {
    a.adj().array()
        += vi.adj().array() * vi.val().array() * digamma(a.val()).array();
  });
}

}  // namespace math
}  // namespace stan
#endif
