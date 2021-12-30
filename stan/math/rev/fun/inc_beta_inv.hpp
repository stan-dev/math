#ifndef STAN_MATH_REV_FUN_INC_BETA_INV_HPP
#define STAN_MATH_REV_FUN_INC_BETA_INV_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/inc_beta_inv.hpp>
#include <stan/math/prim/fun/inc_beta.hpp>
#include <stan/math/prim/fun/exp.hpp>
#include <stan/math/prim/fun/log.hpp>
#include <stan/math/prim/fun/log_diff_exp.hpp>
#include <stan/math/prim/fun/lbeta.hpp>
#include <stan/math/prim/fun/lgamma.hpp>
#include <stan/math/prim/fun/digamma.hpp>
#include <stan/math/prim/fun/F32.hpp>
#include <stan/math/prim/fun/is_any_nan.hpp>

namespace stan {
namespace math {

/**
 * The fused multiply-add function for three variables (C99).
 * This function returns the product of the first two arguments
 * plus the third argument.
 *
 * The partial derivatives are
 *
 * \f$\frac{\partial}{\partial x} (x * y) + z = y\f$, and
 *
 * \f$\frac{\partial}{\partial y} (x * y) + z = x\f$, and
 *
 * \f$\frac{\partial}{\partial z} (x * y) + z = 1\f$.
 *
 * @param x First multiplicand.
 * @param y Second multiplicand.
 * @param z Summand.
 * @return Product of the multiplicands plus the summand, ($a * $b) + $c.
 */
template <typename T1, typename T2, typename T3,
          require_all_stan_scalar_t<T1, T2, T3>* = nullptr,
          require_any_var_t<T1, T2, T3>* = nullptr>
inline var inc_beta_inv(const T1& a, const T2& b, const T3& p) {
  double a_val = value_of(a);
  double b_val = value_of(b);
  double p_val = value_of(p);
  double w = inc_beta_inv(a_val, b_val, p_val);
  return make_callback_var(w, [a, b, p, a_val, b_val, p_val, w](auto& vi) {
    double log_w = log(w);
    double log1m_w = log1m(w);
    double one_m_a = 1 - a_val;
    double one_m_b = 1 - b_val;
    double one_m_w = 1 - w;
    double ap1 = a_val + 1;
    double bp1 = b_val + 1;
    double lbeta_ab = lbeta(a_val, b_val);
    double digamma_apb = digamma(a_val + b_val);

    if(!is_constant_all<T1>::value) {
      double da1 = one_m_b * log1m_w + one_m_a * log_w;
      double da2 = a_val * log_w + 2 * lgamma(a_val)
                   + log(F32(a_val, a_val, one_m_b, ap1, ap1, w))
                   - 2 * lgamma(ap1);
      double da3 = lbeta_ab + log(inc_beta(a_val, b_val, w))
                   + log(log_w - digamma(a_val) + digamma_apb);

      forward_as<var>(a).adj() += vi.adj() * exp(da1 + log_diff_exp(da2, da3));
    }

    if(!is_constant_all<T2>::value) {
      double db1 = (w - 1) * exp(-b_val * log1m_w + one_m_a * log_w);
      double db2 = 2 * lgamma(b_val) 
                   + log(F32(b_val, b_val, one_m_a, bp1, bp1, one_m_w))
                   - 2 * lgamma(bp1)
                   + b_val * log1m_w;

      double db3 = log(inc_beta(b_val, a_val, one_m_w)) + lbeta_ab
                   + log(log1m_w - digamma(b_val) + digamma_apb);

      forward_as<var>(b).adj() += vi.adj() * db1 * exp(log_diff_exp(db2, db3));
    } 

    if(!is_constant_all<T3>::value) {
      forward_as<var>(p).adj() += vi.adj() *
        exp(one_m_b * log1m_w + one_m_a * log_w + lbeta_ab);
    } 
  });
}


}  // namespace math
}  // namespace stan
#endif
