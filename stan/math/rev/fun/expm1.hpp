#ifndef STAN_MATH_REV_FUN_EXPM1_HPP
#define STAN_MATH_REV_FUN_EXPM1_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/prim/fun/expm1.hpp>

namespace stan {
namespace math {

/**
 * The exponentiation of the specified variable minus 1 (C99).
 *
 * The derivative is given by
 *
 * \f$\frac{d}{dx} \exp(a) - 1 = \exp(a)\f$.
 *
 *
   \f[
   \mbox{expm1}(x) =
   \begin{cases}
     e^x-1 & \mbox{if } -\infty\leq x \leq \infty \\[6pt]
     \textrm{NaN} & \mbox{if } x = \textrm{NaN}
   \end{cases}
   \f]

   \f[
   \frac{\partial\, \mbox{expm1}(x)}{\partial x} =
   \begin{cases}
     e^x & \mbox{if } -\infty\leq x\leq \infty \\[6pt]
     \textrm{NaN} & \mbox{if } x = \textrm{NaN}
   \end{cases}
   \f]
 *
 * @param a The variable.
 * @return Two to the power of the specified variable.
 */
inline var expm1(const var& a) {
  return make_callback_var(expm1(a.val()), [a](auto& vi) mutable {
    a.adj() += vi.adj() * (vi.val() + 1.0);
  });
}

template <typename T, require_eigen_t<T>* = nullptr>
inline auto expm1(const var_value<T>& a) {
  return make_callback_var(expm1(a.val()), [a](auto& vi) mutable {
    a.adj().array() += vi.adj().array() * (vi.val().array() + 1.0);
  });
}

}  // namespace math
}  // namespace stan
#endif
