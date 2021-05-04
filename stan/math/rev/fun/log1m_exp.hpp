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
 * @param[in] x Argument.
 * @return Natural logarithm of one minus the exponential of the
 * argument.
 */
inline var log1m_exp(const var& x) {
  return make_callback_var(log1m_exp(x.val()), [x](auto& vi) mutable {
    x.adj() -= vi.adj() / expm1(-x.val());
  });
}

}  // namespace math
}  // namespace stan
#endif
