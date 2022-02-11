#ifndef STAN_MATH_REV_FUN_LOG_HPP
#define STAN_MATH_REV_FUN_LOG_HPP

#include <stan/math/prim/fun/log.hpp>
#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/functor/apply_scalar_unary.hpp>
#include <stan/math/rev/fun/abs.hpp>
#include <stan/math/rev/fun/arg.hpp>
#include <stan/math/rev/fun/atan2.hpp>
#include <stan/math/rev/fun/cos.hpp>
#include <stan/math/rev/fun/is_inf.hpp>
#include <stan/math/rev/fun/is_nan.hpp>
#include <stan/math/rev/fun/norm.hpp>
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
 * @tparam T Arithmetic or a type inheriting from `EigenBase`.
 * @param a Variable whose log is taken.
 * @return Natural log of variable.
 */
template <typename T, require_stan_scalar_or_eigen_t<T>* = nullptr>
inline auto log(const var_value<T>& a) {
  return make_callback_var(log(a.val()), [a](auto& vi) mutable {
    as_array_or_scalar(a.adj())
        += as_array_or_scalar(vi.adj()) / as_array_or_scalar(a.val());
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

}  // namespace math
}  // namespace stan
#endif
