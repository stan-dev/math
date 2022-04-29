#ifndef STAN_MATH_REV_FUN_ABS_HPP
#define STAN_MATH_REV_FUN_ABS_HPP

#include <stan/math/prim/fun/abs.hpp>
#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/fun/fabs.hpp>
#include <stan/math/rev/fun/hypot.hpp>
#include <complex>

namespace stan {
namespace math {

/**
 * Return the absolute value of the variable (std).
 *
 * Delegates to <code>fabs()</code> (see for doc).
 *
   \f[
   \mbox{abs}(x) =
   \begin{cases}
     |x| & \mbox{if } -\infty\leq x\leq \infty \\[6pt]
     \textrm{NaN} & \mbox{if } x = \textrm{NaN}
   \end{cases}
   \f]

   \f[
   \frac{\partial\, \mbox{abs}(x)}{\partial x} =
   \begin{cases}
     -1 & \mbox{if } x < 0 \\
     0 & \mbox{if } x = 0 \\
     1 & \mbox{if } x > 0 \\[6pt]
     \textrm{NaN} & \mbox{if } x = \textrm{NaN}
   \end{cases}
   \f]
 *
 * @tparam T A floating point type or an Eigen type with floating point scalar.
 * @param a Variable input.
 * @return Absolute value of variable.
 */
template <typename T>
inline auto abs(const var_value<T>& a) {
  return fabs(a);
}

/**
 * Return the absolute value of the complex argument.
 *
 * @param[in] z argument
 * @return absolute value of the argument
 */
inline var abs(const std::complex<var>& z) { return internal::complex_abs(z); }

}  // namespace math
}  // namespace stan
#endif
