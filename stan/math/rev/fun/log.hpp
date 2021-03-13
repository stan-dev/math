#ifndef STAN_MATH_REV_FUN_LOG_HPP
#define STAN_MATH_REV_FUN_LOG_HPP

#include <stan/math/prim/fun/log.hpp>
#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/fun/abs.hpp>
#include <stan/math/rev/fun/arg.hpp>
#include <stan/math/rev/fun/atan2.hpp>
#include <stan/math/rev/fun/cos.hpp>
#include <stan/math/rev/fun/is_inf.hpp>
#include <stan/math/rev/fun/is_nan.hpp>
#include <stan/math/rev/fun/sqrt.hpp>
#include <cmath>

namespace stan {
namespace math {

/**
 * Return the natural log of the specified variable (cmath).
 *
 * The derivative is defined by
 *
 * \f$\frac{d}{dx} \log x = \frac{1}{x}\f$.
 *
   \f[
   \mbox{log}(x) =
   \begin{cases}
     \textrm{NaN} & \mbox{if } x < 0\\
     \ln(x) & \mbox{if } x \geq 0 \\[6pt]
     \textrm{NaN} & \mbox{if } x = \textrm{NaN}
   \end{cases}
   \f]

   \f[
   \frac{\partial\, \mbox{log}(x)}{\partial x} =
   \begin{cases}
     \textrm{NaN} & \mbox{if } x < 0\\
     \frac{1}{x} & \mbox{if } x\geq 0 \\[6pt]
     \textrm{NaN} & \mbox{if } x = \textrm{NaN}
   \end{cases}
   \f]
 *
 * @param a Variable whose log is taken.
 * @return Natural log of variable.
 */
inline var log(const var& a) {
  return make_callback_var(std::log(a.val()), [a](auto& vi) mutable {
    a.adj() += vi.adj() / a.val();
  });
}

/**
 * Return the natural logarithm (base e) of the specified complex argument.
 *
 * @param z complex argument
 * @return natural logarithm of argument
 */
inline std::complex<var> log(const std::complex<var>& z) {
  return internal::complex_log(z);
}

/**
 * Return the natural log of the elements of x
 *
 * @tparam T type of x
 * @param x argument
 * @return elementwise natural log of x
 */
template <typename T, require_var_matrix_t<T>* = nullptr>
inline auto log(const T& x) {
  return make_callback_var(
      x.val().array().log().matrix(), [x](const auto& vi) mutable {
        x.adj() += (vi.adj().array() / x.val().array()).matrix();
      });
}

}  // namespace math
}  // namespace stan
#endif
