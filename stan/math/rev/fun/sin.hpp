#ifndef STAN_MATH_REV_FUN_SIN_HPP
#define STAN_MATH_REV_FUN_SIN_HPP

#include <stan/math/prim/fun/cos.hpp>
#include <stan/math/prim/fun/isfinite.hpp>
#include <stan/math/prim/fun/isinf.hpp>
#include <stan/math/prim/fun/sin.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/fun/is_inf.hpp>
#include <stan/math/rev/fun/cosh.hpp>
#include <stan/math/rev/fun/sinh.hpp>
#include <cmath>
#include <complex>

namespace stan {
namespace math {

/**
 * Return the sine of a radian-scaled variable (cmath).
 *
 * The derivative is defined by
 *
 * \f$\frac{d}{dx} \sin x = \cos x\f$.
 *
 *
   \f[
   \mbox{sin}(x) =
   \begin{cases}
     \sin(x) & \mbox{if } -\infty\leq x \leq \infty \\[6pt]
     \textrm{NaN} & \mbox{if } x = \textrm{NaN}
   \end{cases}
   \f]

   \f[
   \frac{\partial\, \mbox{sin}(x)}{\partial x} =
   \begin{cases}
     \cos(x) & \mbox{if } -\infty\leq x\leq \infty \\[6pt]
     \textrm{NaN} & \mbox{if } x = \textrm{NaN}
   \end{cases}
   \f]
 *
 * @param a Variable for radians of angle.
 * @return Sine of variable.
 */
inline var sin(const var& a) {
  return make_callback_var(std::sin(a.val()), [a](const auto& vi) mutable {
    a.adj() += vi.adj() * std::cos(a.val());
  });
}

/**
 * Return the sine of a radian-scaled variable (cmath).
 *
 * @tparam Varmat a `var_value` with inner Eigen type
 * @param a Variable for radians of angle.
 * @return Sine of variable.
 */
template <typename VarMat, require_var_matrix_t<VarMat>* = nullptr>
inline auto sin(const VarMat& a) {
  return make_callback_var(
      a.val().array().sin().matrix(), [a](const auto& vi) mutable {
        a.adj() += vi.adj().cwiseProduct(a.val().array().cos().matrix());
      });
}

/**
 * Return the sine of the complex argument.
 *
 * @param[in] z argument
 * @return sine of the argument
 */
inline std::complex<var> sin(const std::complex<var>& z) {
  return stan::math::internal::complex_sin(z);
}

}  // namespace math
}  // namespace stan
#endif
