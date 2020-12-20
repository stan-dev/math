#ifndef STAN_MATH_REV_FUN_SINH_HPP
#define STAN_MATH_REV_FUN_SINH_HPP

#include <stan/math/prim/core.hpp>
#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/isinf.hpp>
#include <stan/math/prim/fun/sinh.hpp>
#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/fun/exp.hpp>
#include <cmath>
#include <complex>

namespace stan {
namespace math {

/**
 * Return the hyperbolic sine of the specified variable (cmath).
 *
 * The derivative is defined by
 *
 * \f$\frac{d}{dx} \sinh x = \cosh x\f$.
 *
 *
   \f[
   \mbox{sinh}(x) =
   \begin{cases}
     \sinh(x) & \mbox{if } -\infty\leq x \leq \infty \\[6pt]
     \textrm{NaN} & \mbox{if } x = \textrm{NaN}
   \end{cases}
   \f]

   \f[
   \frac{\partial\, \mbox{sinh}(x)}{\partial x} =
   \begin{cases}
     \cosh(x) & \mbox{if } -\infty\leq x\leq \infty \\[6pt]
     \textrm{NaN} & \mbox{if } x = \textrm{NaN}
   \end{cases}
   \f]
 *
 * @param a Variable.
 * @return Hyperbolic sine of variable.
 */
inline var sinh(const var& a) {
  return make_callback_var(std::sinh(a.val()), [a](const auto& vi) mutable {
    a.adj() += vi.adj() * std::cosh(a.val());
  });
}

/**
 * Return the hyperbolic of a radian-scaled variable (cmath).
 *
 *
 * @tparam Varmat a `var_value` with inner Eigen type
 * @param a Variable for radians of angle.
 * @return Hyperbolid Sine of variable.
 */
template <typename VarMat, require_var_matrix_t<VarMat>* = nullptr>
inline auto sinh(const VarMat& a) {
  return make_callback_var(
      a.val().array().sinh().matrix(), [a](const auto& vi) mutable {
        a.adj() += vi.adj().cwiseProduct(a.val().array().cosh().matrix());
      });
}

/**
 * Return the hyperbolic sine of the complex argument.
 *
 * @param[in] z argument
 * @return hyperbolic sine of the argument
 */
inline std::complex<var> sinh(const std::complex<var>& z) {
  return stan::math::internal::complex_sinh(z);
}

}  // namespace math
}  // namespace stan
#endif
