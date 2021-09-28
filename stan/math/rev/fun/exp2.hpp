#ifndef STAN_MATH_REV_FUN_EXP2_HPP
#define STAN_MATH_REV_FUN_EXP2_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/prim/fun/exp2.hpp>
#include <stan/math/prim/fun/constants.hpp>
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
inline var exp2(const var& a) {
  return make_callback_var(std::exp2(a.val()), [a](auto& vi) mutable {
    a.adj() += vi.adj() * vi.val() * LOG_TWO;
  });
}

template <typename T, require_eigen_t<T>* = nullptr>
inline auto exp2(const var_value<T>& a) {
  return make_callback_var(exp2(a.val()), [a](auto& vi) mutable {
    a.adj().array() += vi.adj().array() * vi.val().array() * LOG_TWO;
  });
}

}  // namespace math
}  // namespace stan
#endif
