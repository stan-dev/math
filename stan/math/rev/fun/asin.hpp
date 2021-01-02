#ifndef STAN_MATH_REV_FUN_ASIN_HPP
#define STAN_MATH_REV_FUN_ASIN_HPP

#include <stan/math/prim/fun/asin.hpp>
#include <stan/math/prim/fun/abs.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/fun/abs.hpp>
#include <stan/math/rev/fun/asinh.hpp>
#include <stan/math/rev/fun/value_of_rec.hpp>
#include <cmath>
#include <complex>

namespace stan {
namespace math {

/**
 * Return the principal value of the arc sine, in radians, of the
 * specified variable (cmath).
 *
 * The derivative is defined by
 *
 * \f$\frac{d}{dx} \arcsin x = \frac{1}{\sqrt{1 - x^2}}\f$.
 *
 *
   \f[
   \mbox{asin}(x) =
   \begin{cases}
     \textrm{NaN} & \mbox{if } x < -1\\
     \arcsin(x) & \mbox{if } -1\leq x\leq 1 \\
     \textrm{NaN} & \mbox{if } x > 1\\[6pt]
     \textrm{NaN} & \mbox{if } x = \textrm{NaN}
   \end{cases}
   \f]

   \f[
   \frac{\partial\, \mbox{asin}(x)}{\partial x} =
   \begin{cases}
     \textrm{NaN} & \mbox{if } x < -1\\
     \frac{\partial\, \arcsin(x)}{\partial x} & \mbox{if } -1\leq x\leq 1 \\
     \textrm{NaN} & \mbox{if } x < -1\\[6pt]
     \textrm{NaN} & \mbox{if } x = \textrm{NaN}
   \end{cases}
   \f]

   \f[
   \frac{\partial \, \arcsin(x)}{\partial x} = \frac{1}{\sqrt{1-x^2}}
   \f]
 *
 * @param x Variable in range [-1, 1].
 * @return Arc sine of variable, in radians.
 */
inline var asin(const var& x) {
  return make_callback_var(std::asin(x.val()), [x](const auto& vi) mutable {
    x.adj() += vi.adj() / std::sqrt(1.0 - (x.val() * x.val()));
  });
}

/**
 * Return the principal value of the arc sine, in radians, of the
 * specified variable (cmath).
 *
 * @tparam Varmat a `var_value` with inner Eigen type
 * @param x Variable with cells in range [-1, 1].
 * @return Arc sine of variable, in radians.
 */
template <typename VarMat, require_var_matrix_t<VarMat>* = nullptr>
inline auto asin(const VarMat& x) {
  return make_callback_var(
      x.val().array().asin().matrix(), [x](const auto& vi) mutable {
        x.adj().array()
            += vi.adj().array() / (1.0 - (x.val().array().square())).sqrt();
      });
}

/**
 * Return the arc sine of the complex argument.
 *
 * @param[in] z argument
 * @return arc sine of the argument
 */
inline std::complex<var> asin(const std::complex<var>& z) {
  return stan::math::internal::complex_asin(z);
}

}  // namespace math
}  // namespace stan
#endif
