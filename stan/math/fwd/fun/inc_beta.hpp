#ifndef STAN_MATH_FWD_FUN_INC_BETA_HPP
#define STAN_MATH_FWD_FUN_INC_BETA_HPP

#include <stan/math/fwd/core.hpp>
#include <stan/math/fwd/meta.hpp>
#include <stan/math/fwd/fun/digamma.hpp>
#include <stan/math/prim/fun/grad_reg_inc_beta.hpp>
#include <stan/math/prim/fun/inc_beta.hpp>
#include <stan/math/prim/fun/inc_beta_ddz.hpp>
#include <stan/math/prim/fun/value_of.hpp>

namespace stan {
namespace math {

/**
 * The regularized incomplete beta function I_x(a, b) for <code>fvar</code>
 * arguments. Only the partials of the <code>fvar</code> arguments are
 * computed. Each argument is an <code>fvar</code> or arithmetic; mixing with
 * <code>var</code> is not supported.
 *
 * @tparam T1 type of a
 * @tparam T2 type of b
 * @tparam T3 type of x
 * @param a first shape, positive
 * @param b second shape, positive
 * @param x upper limit, in [0, 1]
 * @return I_x(a, b)
 */
template <typename T1, typename T2, typename T3,
          require_all_stan_scalar_t<T1, T2, T3>* = nullptr,
          require_all_not_var_t<T1, T2, T3>* = nullptr,
          require_any_fvar_t<T1, T2, T3>* = nullptr>
inline fvar<partials_return_t<T1, T2, T3>> inc_beta(const T1& a, const T2& b,
                                                    const T3& x) {
  using T_return = partials_return_t<T1, T2, T3>;
  const auto a_val = value_of(a);
  const auto b_val = value_of(b);
  const auto x_val = value_of(x);
  T_return d(0);
  if constexpr (is_autodiff_v<T1> && is_autodiff_v<T2>) {
    const auto [d_a, d_b] = grad_reg_inc_beta<GRAD_RIB::BOTH>(
        a_val, b_val, x_val, digamma(a_val), digamma(b_val),
        digamma(a_val + b_val), 0.0);
    d += a.d_ * d_a + b.d_ * d_b;
  } else if constexpr (is_autodiff_v<T1>) {
    d += a.d_
         * grad_reg_inc_beta<GRAD_RIB::G1>(a_val, b_val, x_val, digamma(a_val),
                                           0.0, digamma(a_val + b_val), 0.0);
  } else if constexpr (is_autodiff_v<T2>) {
    d += b.d_
         * grad_reg_inc_beta<GRAD_RIB::G2>(a_val, b_val, x_val, 0.0,
                                           digamma(b_val),
                                           digamma(a_val + b_val), 0.0);
  }
  if constexpr (is_autodiff_v<T3>) {
    d += x.d_ * inc_beta_ddz(a_val, b_val, x_val);
  }
  return fvar<T_return>(inc_beta(a_val, b_val, x_val), d);
}

}  // namespace math
}  // namespace stan

#endif
