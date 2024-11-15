#ifndef STAN_MATH_REV_FUN_COS_HPP
#define STAN_MATH_REV_FUN_COS_HPP

#include <stan/math/rev/core.hpp>
#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/fun/abs.hpp>
#include <stan/math/rev/fun/cosh.hpp>
#include <stan/math/rev/fun/sinh.hpp>
#include <stan/math/prim/fun/cos.hpp>
#include <cmath>
#include <complex>

namespace stan {
namespace math {

/**
 * Return the cosine of a radian-scaled variable (cmath).
 *
 * The derivative is defined by
 *
 * \f$\frac{d}{dx} \cos x = - \sin x\f$.
 *
 *
   \f[
   \mbox{cos}(x) =
   \begin{cases}
     \cos(x) & \mbox{if } -\infty\leq x \leq \infty \\[6pt]
     \textrm{NaN} & \mbox{if } x = \textrm{NaN}
   \end{cases}
   \f]

   \f[
   \frac{\partial\, \mbox{cos}(x)}{\partial x} =
   \begin{cases}
     -\sin(x) & \mbox{if } -\infty\leq x\leq \infty \\[6pt]
     \textrm{NaN} & \mbox{if } x = \textrm{NaN}
   \end{cases}
   \f]
 *
 * @param a Variable for radians of angle.
 * @return Cosine of variable.
 */
inline var cos(var a) {
  return make_callback_var(std::cos(a.val()), [a](const auto& vi) mutable {
    a.adj() -= vi.adj() * std::sin(a.val());
  });
}

/**
 * Return the cosine of a radian-scaled variable (cmath).
 *
 *
 * @tparam Varmat a `var_value` with inner Eigen type
 * @param a Variable for radians of angle.
 * @return Cosine of variable.
 */
template <typename VarMat, require_var_matrix_t<VarMat>* = nullptr>
inline auto cos(const VarMat& a) {
  return make_callback_var(
      a.val().array().cos().matrix(), [a](const auto& vi) mutable {
        a.adj() -= vi.adj().cwiseProduct(a.val().array().sin().matrix());
      });
}
/**
 * Return the cosine of the complex argument.
 *
 * @param[in] z argument
 * @return cosine of the argument
 */
inline std::complex<var> cos(const std::complex<var>& z) {
  return stan::math::internal::complex_cos(z);
}

}  // namespace math
}  // namespace stan
#endif
