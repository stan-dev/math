#ifndef STAN_MATH_REV_FUN_EXP2_HPP
#define STAN_MATH_REV_FUN_EXP2_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/exp2.hpp>
#include <cmath>

namespace stan {
namespace math {

/**
 * Exponentiation base 2 function for variables (C99).
 *
 * The derivative is
 *
 * \f$\frac{d}{dx} 2^x = (\log 2) 2^x\f$.
 *
   \f[
   \mbox{exp2}(x) =
   \begin{cases}
     2^x & \mbox{if } -\infty\leq x \leq \infty \\[6pt]
     \textrm{NaN} & \mbox{if } x = \textrm{NaN}
   \end{cases}
   \f]

   \f[
   \frac{\partial\, \mbox{exp2}(x)}{\partial x} =
   \begin{cases}
     2^x\ln2 & \mbox{if } -\infty\leq x\leq \infty \\[6pt]
     \textrm{NaN} & \mbox{if } x = \textrm{NaN}
   \end{cases}
   \f]
 *
 * @param a The variable.
 * @return Two to the power of the specified variable.
 */
template <typename T, require_var_t<T>* = nullptr>
inline auto exp2(T&& a) {
  return make_callback_var(exp2(a.val()), [a](auto& vi) mutable {
    as_array_or_scalar(a.adj()) += as_array_or_scalar(vi.adj())
                                   * as_array_or_scalar(vi.val()) * LOG_TWO;
  });
}

}  // namespace math
}  // namespace stan
#endif
