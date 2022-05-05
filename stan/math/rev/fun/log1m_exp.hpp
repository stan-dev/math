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
inline auto log1m_exp(const var& x) {
  return make_callback_var(log1m_exp(x.val()), [x](auto& vi) mutable {
    x.adj() -= vi.adj() / expm1(-x.val());
  });
}

template <typename T, require_rev_matrix_t<T>* = nullptr>
inline auto log1m_exp(const T& x) {
  auto x_arena = to_arena(x);
  return make_callback_rev_matrix<T>(log1m_exp(x_arena.val()), [x_arena](auto&& vi) mutable {
    x_arena.adj().array()
        -= vi.adj().array() / expm1(-x_arena.val()).array();
  });
}

}  // namespace math
}  // namespace stan
#endif
