#ifndef STAN_MATH_REV_FUN_TAN_HPP
#define STAN_MATH_REV_FUN_TAN_HPP

#include <stan/math/prim/fun/isinf.hpp>
#include <stan/math/prim/fun/tan.hpp>
#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/fun/is_inf.hpp>
#include <stan/math/rev/fun/sinh.hpp>
#include <stan/math/rev/fun/tanh.hpp>
#include <cmath>
#include <complex>

namespace stan {
namespace math {

/**
 * Return the tangent of a radian-scaled variable (cmath).
 *
 * The derivative is defined by
 *
 * \f$\frac{d}{dx} \tan x = \sec^2 x\f$.
 *
 *
   \f[
   \mbox{tan}(x) =
   \begin{cases}
     \tan(x) & \mbox{if } -\infty\leq x \leq \infty \\[6pt]
     \textrm{NaN} & \mbox{if } x = \textrm{NaN}
   \end{cases}
   \f]

   \f[
   \frac{\partial\, \mbox{tan}(x)}{\partial x} =
   \begin{cases}
     \sec^2(x) & \mbox{if } -\infty\leq x\leq \infty \\[6pt]
     \textrm{NaN} & \mbox{if } x = \textrm{NaN}
   \end{cases}
   \f]
 *
 * @param a Variable for radians of angle.
 * @return Tangent of variable.
 */
inline var tan(const var& a) {
  return make_callback_var(std::tan(a.val()), [a](const auto& vi) mutable {
    a.adj() += vi.adj() * (1.0 + vi.val() * vi.val());
  });
}

/**
 * Return the tangent of a radian-scaled variable (cmath).
 *
 *
 * @tparam Varmat a `var_value` with inner Eigen type
 * @param a Variable for radians of angle.
 * @return Tangent of variable.
 */
template <typename VarMat, require_var_matrix_t<VarMat>* = nullptr>
inline auto tan(const VarMat& a) {
  return make_callback_var(a.val().array().tan().matrix(),
                           [a](const auto& vi) mutable {
                             a.adj() += vi.adj().cwiseProduct(
                                 (1.0 + vi.val().array().square()).matrix());
                           });
}

/**
 * Return the tangent of the complex argument.
 *
 * @param[in] z argument
 * @return tangent of the argument
 */
inline std::complex<var> tan(const std::complex<var>& z) {
  return stan::math::internal::complex_tan(z);
}

}  // namespace math
}  // namespace stan
#endif
