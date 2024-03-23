#ifndef STAN_MATH_REV_FUN_SQRT_HPP
#define STAN_MATH_REV_FUN_SQRT_HPP

#include <stan/math/prim/fun/sqrt.hpp>
#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/fun/atan2.hpp>
#include <stan/math/rev/fun/cos.hpp>
#include <stan/math/rev/fun/hypot.hpp>
#include <cmath>
#include <complex>

namespace stan {
namespace math {

/**
 * Return the square root of the specified variable (cmath).
 *
 * The derivative is defined by
 *
 * \f$\frac{d}{dx} \sqrt{x} = \frac{1}{2 \sqrt{x}}\f$.
 *
   \f[
   \mbox{sqrt}(x) =
   \begin{cases}
     \textrm{NaN} & x < 0 \\
     \sqrt{x} & \mbox{if } x\geq 0\\[6pt]
     \textrm{NaN} & \mbox{if } x = \textrm{NaN}
   \end{cases}
   \f]

   \f[
   \frac{\partial\, \mbox{sqrt}(x)}{\partial x} =
   \begin{cases}
     \textrm{NaN} & x < 0 \\
     \frac{1}{2\sqrt{x}} & x\geq 0\\[6pt]
     \textrm{NaN} & \mbox{if } x = \textrm{NaN}
   \end{cases}
   \f]
 *
 * @param a Variable whose square root is taken.
 * @return Square root of variable.
 */
inline var sqrt(const var& a) {
  return make_callback_var(std::sqrt(a.val()), [a](auto& vi) mutable {
    if (vi.val() != 0.0) {
      a.adj() += vi.adj() / (2.0 * vi.val());
    }
  });
}

/**
 * Return elementwise square root of vector
 *
 * @tparam T a `var_value` with inner Eigen type
 * @param a input
 * @return elementwise square root of vector
 */
template <typename T, require_var_matrix_t<T>* = nullptr>
inline auto sqrt(const T& a) {
  return make_callback_var(
      a.val().array().sqrt().matrix(), [a](auto& vi) mutable {
        a.adj().array()
            += (vi.val_op().array() == 0.0)
                   .select(0.0, vi.adj().array() / (2.0 * vi.val_op().array()));
      });
}

/**
 * Return the square root of the complex argument.
 *
 * @param[in] z argument
 * @return square root of the argument
 */
inline std::complex<var> sqrt(const std::complex<var>& z) {
  return internal::complex_sqrt(z);
}

}  // namespace math
}  // namespace stan
#endif
