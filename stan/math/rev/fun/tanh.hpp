#ifndef STAN_MATH_REV_FUN_TANH_HPP
#define STAN_MATH_REV_FUN_TANH_HPP

#include <stan/math/prim/fun/tanh.hpp>
#include <stan/math/prim/fun/cosh.hpp>
#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/fun/exp.hpp>
#include <cmath>
#include <complex>

namespace stan {
namespace math {

/**
 * Return the hyperbolic tangent of the specified variable (cmath).
 *
 * The derivative is defined by
 *
 * \f$\frac{d}{dx} \tanh x = \frac{1}{\cosh^2 x}\f$.
 *
 *
   \f[
   \mbox{tanh}(x) =
   \begin{cases}
     \tanh(x) & \mbox{if } -\infty\leq x \leq \infty \\[6pt]
     \textrm{NaN} & \mbox{if } x = \textrm{NaN}
   \end{cases}
   \f]

   \f[
   \frac{\partial\, \mbox{tanh}(x)}{\partial x} =
   \begin{cases}
     \mbox{sech}^2(x) & \mbox{if } -\infty\leq x\leq \infty \\[6pt]
     \textrm{NaN} & \mbox{if } x = \textrm{NaN}
   \end{cases}
   \f]
 *
 * @param a Variable.
 * @return Hyperbolic tangent of variable.
 */
inline var tanh(const var& a) {
  return make_callback_var(std::tanh(a.val()), [a](const auto& vi) mutable {
    const auto a_cosh = std::cosh(a.val());
    a.adj() += vi.adj_ / (a_cosh * a_cosh);
  });
}

/**
 * Return the hyperbolic tangent of elements of a
 *
 * @tparam T type of a
 * @param a argument
 * @return elementwise hyperbolic tangent of a
 */
template <typename VarMat, require_var_matrix_t<VarMat>* = nullptr>
inline auto tanh(const VarMat& a) {
  return make_callback_var(
      a.val().array().tanh().matrix(), [a](const auto& vi) mutable {
        a.adj().array() += vi.adj_.array() / (a.val().array().cosh().square());
      });
}

/**
 * Return the hyperbolic tangent of the complex argument.
 *
 * @param[in] z argument
 * @return hyperbolic tangent of the argument
 */
inline std::complex<var> tanh(const std::complex<var>& z) {
  return stan::math::internal::complex_tanh(z);
}

}  // namespace math
}  // namespace stan
#endif
