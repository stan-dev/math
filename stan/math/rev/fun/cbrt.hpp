#ifndef STAN_MATH_REV_FUN_CBRT_HPP
#define STAN_MATH_REV_FUN_CBRT_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/prim/fun/cbrt.hpp>

namespace stan {
namespace math {

/**
 * Returns the cube root of the specified variable (C99).
 *
 * The derivative is
 *
 * \f$\frac{d}{dx} x^{1/3} = \frac{1}{3 x^{2/3}}\f$.
 *
   \f[
   \mbox{cbrt}(x) =
   \begin{cases}
     \sqrt[3]{x} & \mbox{if } -\infty\leq x \leq \infty \\[6pt]
     \textrm{NaN} & \mbox{if } x = \textrm{NaN}
   \end{cases}
   \f]

   \f[
   \frac{\partial\, \mbox{cbrt}(x)}{\partial x} =
   \begin{cases}
     \frac{1}{3x^{2/3}} & \mbox{if } -\infty\leq x\leq \infty \\[6pt]
     \textrm{NaN} & \mbox{if } x = \textrm{NaN}
   \end{cases}
   \f]
 *
 * @param a Specified variable.
 * @return Cube root of the variable.
 */
inline var cbrt(const var& a) {
  return make_callback_var(cbrt(a.val()), [a](const auto& vi) mutable {
     a.adj() += vi.adj_ / (3.0 * vi.val_ * vi.val_);
  });
}

template <typename VarMat, require_var_matrix_t<VarMat>* = nullptr>
inline auto cbrt(const VarMat& a) {
  return make_callback_var(a.val().unaryExpr([](const auto x) {
     return cbrt(x);
   }),
   [a](const auto& vi) mutable {
     a.adj().array() += vi.adj_.array() / (3.0 * vi.val_.array().square());
  });
}


}  // namespace math
}  // namespace stan
#endif
