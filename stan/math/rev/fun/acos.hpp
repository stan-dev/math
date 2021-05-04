#ifndef STAN_MATH_REV_FUN_ACOS_HPP
#define STAN_MATH_REV_FUN_ACOS_HPP

#include <stan/math/prim/fun/acos.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/fun/abs.hpp>
#include <stan/math/rev/fun/arg.hpp>
#include <stan/math/rev/fun/asin.hpp>
#include <stan/math/rev/fun/is_inf.hpp>
#include <stan/math/rev/fun/is_nan.hpp>
#include <stan/math/rev/fun/polar.hpp>
#include <cmath>
#include <complex>

namespace stan {
namespace math {

/**
 * Return the principal value of the arc cosine of a variable,
 * in radians (cmath).
 *
 * The derivative is defined by
 *
 * \f$\frac{d}{dx} \arccos x = \frac{-1}{\sqrt{1 - x^2}}\f$.
 *
 *
   \f[
   \mbox{acos}(x) =
   \begin{cases}
     \textrm{NaN} & \mbox{if } x < -1\\
     \arccos(x) & \mbox{if } -1\leq x\leq 1 \\
     \textrm{NaN} & \mbox{if } x > 1\\[6pt]
     \textrm{NaN} & \mbox{if } x = \textrm{NaN}
   \end{cases}
   \f]

   \f[
   \frac{\partial\, \mbox{acos}(x)}{\partial x} =
   \begin{cases}
     \textrm{NaN} & \mbox{if } x < -1\\
     \frac{\partial\, \arccos(x)}{\partial x} & \mbox{if } -1\leq x\leq 1 \\
     \textrm{NaN} & \mbox{if } x < -1\\[6pt]
     \textrm{NaN} & \mbox{if } x = \textrm{NaN}
   \end{cases}
   \f]

   \f[
   \frac{\partial \, \arccos(x)}{\partial x} = -\frac{1}{\sqrt{1-x^2}}
   \f]
 *
 * @param x argument
 * @return Arc cosine of variable, in radians.
 */
inline var acos(const var& x) {
  return make_callback_var(std::acos(x.val()), [x](const auto& vi) mutable {
    x.adj() -= vi.adj() / std::sqrt(1.0 - (x.val() * x.val()));
  });
}

/**
 * Return the principal value of the arc cosine of a variable,
 * in radians (cmath).
 *
 * @param x a `var_value` with inner Eigen type
 * @return Arc cosine of variable, in radians.
 */
template <typename VarMat, require_var_matrix_t<VarMat>* = nullptr>
inline auto acos(const VarMat& x) {
  return make_callback_var(
      x.val().array().acos().matrix(), [x](const auto& vi) mutable {
        x.adj().array()
            -= vi.adj().array() / (1.0 - (x.val().array().square())).sqrt();
      });
}

/**
 * Return the arc cosine of the complex argument.
 *
 * @param x argument
 * @return arc cosine of the argument
 */
inline std::complex<var> acos(const std::complex<var>& x) {
  return stan::math::internal::complex_acos(x);
}

}  // namespace math
}  // namespace stan
#endif
