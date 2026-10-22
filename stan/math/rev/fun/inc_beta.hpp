#ifndef STAN_MATH_REV_FUN_INC_BETA_HPP
#define STAN_MATH_REV_FUN_INC_BETA_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/prim/fun/digamma.hpp>
#include <stan/math/prim/fun/grad_reg_inc_beta.hpp>
#include <stan/math/prim/fun/inc_beta.hpp>
#include <stan/math/prim/fun/inc_beta_ddz.hpp>
#include <stan/math/prim/fun/value_of.hpp>

namespace stan {
namespace math {

/**
 * The regularized incomplete beta function I_z(a, b) for <code>var</code>
 * arguments. Only the partials of the <code>var</code> arguments are
 * computed in the reverse pass. Each argument is a <code>var</code> or
 * arithmetic. Mixing with <code>fvar</code> is not supported.
 *
 * @tparam T1 type of a
 * @tparam T2 type of b
 * @tparam T3 type of z
 * @param a first shape, positive
 * @param b second shape, positive
 * @param z upper limit, in [0, 1]
 * @return I_z(a, b)
 */
template <typename T1, typename T2, typename T3,
          require_all_stan_scalar_t<T1, T2, T3>* = nullptr,
          require_all_var_or_arithmetic_t<T1, T2, T3>* = nullptr,
          require_any_var_t<T1, T2, T3>* = nullptr>
inline var inc_beta(const T1& a, const T2& b, const T3& z) {
  return make_callback_var(
      inc_beta(value_of(a), value_of(b), value_of(z)), [a, b, z](auto& vi) {
        const double a_val = value_of(a);
        const double b_val = value_of(b);
        const double z_val = value_of(z);
        if constexpr (is_autodiff_v<T1> && is_autodiff_v<T2>) {
          const auto [d_a, d_b] = grad_reg_inc_beta<GRAD_RIB::BOTH>(
              a_val, b_val, z_val, digamma(a_val), digamma(b_val),
              digamma(a_val + b_val), 0.0);
          a.adj() += vi.adj() * d_a;
          b.adj() += vi.adj() * d_b;
        } else if constexpr (is_autodiff_v<T1>) {
          a.adj() += vi.adj()
                     * grad_reg_inc_beta<GRAD_RIB::G1>(
                         a_val, b_val, z_val, digamma(a_val), 0.0,
                         digamma(a_val + b_val), 0.0);
        } else if constexpr (is_autodiff_v<T2>) {
          b.adj() += vi.adj()
                     * grad_reg_inc_beta<GRAD_RIB::G2>(
                         a_val, b_val, z_val, 0.0, digamma(b_val),
                         digamma(a_val + b_val), 0.0);
        }
        if constexpr (is_autodiff_v<T3>) {
          z.adj() += vi.adj() * inc_beta_ddz(a_val, b_val, z_val);
        }
      });
}

}  // namespace math
}  // namespace stan
#endif
