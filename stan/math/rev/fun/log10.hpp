#ifndef STAN_MATH_REV_FUN_LOG10_HPP
#define STAN_MATH_REV_FUN_LOG10_HPP

#include <stan/math/prim/core.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/log10.hpp>
#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/fun/atan2.hpp>
#include <stan/math/rev/fun/hypot.hpp>
#include <stan/math/rev/fun/log.hpp>
#include <cmath>
#include <complex>

namespace stan {
namespace math {

/**
 * Return the base 10 log of the specified variable (cmath).
 *
 * The derivative is defined by
 *
 * \f$\frac{d}{dx} \log_{10} x = \frac{1}{x \log 10}\f$.
 *
 *
   \f[
   \mbox{log10}(x) =
   \begin{cases}
     \textrm{NaN} & \mbox{if } x < 0\\
     \log_{10}(x) & \mbox{if } x \geq 0 \\[6pt]
     \textrm{NaN} & \mbox{if } x = \textrm{NaN}
   \end{cases}
   \f]

   \f[
   \frac{\partial\, \mbox{log10}(x)}{\partial x} =
   \begin{cases}
     \textrm{NaN} & \mbox{if } x < 0\\
     \frac{1}{x \ln10} & \mbox{if } x\geq 0 \\[6pt]
     \textrm{NaN} & \mbox{if } x = \textrm{NaN}
   \end{cases}
   \f]
 *
 * @tparam T Arithmetic or a type inheriting from `EigenBase`.
 * @param a Variable whose log is taken.
 * @return Base 10 log of variable.
 */
inline auto log10(const var a) {
  return make_callback_var(log10(a.val()), [a](auto& vi) mutable {
    a.adj() += vi.adj() / (LOG_TEN * a.val());
  });
}

template <typename T, require_rev_matrix_t<T>* = nullptr>
inline auto log10(const T& a) {
  auto a_arena = to_arena(a);
  return make_callback_rev_matrix<T>(log10(a_arena.val()), [a_arena](auto&& vi) mutable {
    a_arena.adj().array() += vi.adj().array() / (LOG_TEN * a_arena.val().array());
  });
}

/**
 * Return the base 10 logarithm of the specified complex number.
 *
 * @param z complex argument
 * @return base 10 log of argument
 */
inline std::complex<var> log10(const std::complex<var>& z) {
  return internal::complex_log10(z);
}

}  // namespace math
}  // namespace stan
#endif
