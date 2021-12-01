#ifndef STAN_MATH_REV_FUN_LOG1M_EXP_HPP
#define STAN_MATH_REV_FUN_LOG1M_EXP_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/prim/fun/expm1.hpp>
#include <stan/math/prim/fun/log1m_exp.hpp>

namespace stan {
namespace math {

/**
 * Return the log of 1 minus the exponential of the specified
 * variable.
 *
 * <p>The derivative of <code>log(1 - exp(x))</code> with respect
 * to <code>x</code> is <code>-1 / expm1(-x)</code>.
 *
 * @tparam T Arithmetic or a type inheriting from `EigenBase`.
 * @param[in] x Argument.
 * @return Natural logarithm of one minus the exponential of the
 * argument.
 */
template <typename T, require_stan_scalar_or_eigen_t<T>* = nullptr>
inline auto log1m_exp(const var_value<T>& x) {
  return make_callback_var(log1m_exp(x.val()), [x](auto& vi) mutable {
    as_array_or_scalar(x.adj())
        -= as_array_or_scalar(vi.adj()) / as_array_or_scalar(expm1(-x.val()));
  });
}

}  // namespace math
}  // namespace stan
#endif
