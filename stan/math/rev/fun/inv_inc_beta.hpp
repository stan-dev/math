#ifndef STAN_MATH_REV_FUN_INV_INC_BETA_HPP
#define STAN_MATH_REV_FUN_INV_INC_BETA_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/fun/fabs.hpp>
#include <stan/math/rev/fun/inc_beta.hpp>
#include <stan/math/rev/fun/digamma.hpp>
#include <stan/math/rev/fun/exp.hpp>
#include <stan/math/rev/fun/hypergeometric_pFq.hpp>
#include <stan/math/rev/fun/log.hpp>
#include <stan/math/rev/fun/log_diff_exp.hpp>
#include <stan/math/rev/fun/lbeta.hpp>
#include <stan/math/rev/fun/lgamma.hpp>
#include <stan/math/rev/fun/sum.hpp>
#include <stan/math/prim/fun/is_any_nan.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/hypergeometric_3F2.hpp>
#include <stan/math/prim/fun/inv_inc_beta.hpp>

namespace stan {
namespace math {

/**
 * The inverse of the normalized incomplete beta function of a, b, with
 * probability p.
 *
 * Used to compute the inverse cumulative density function for the beta
 * distribution.
 *
   \f[
   \frac{\partial }{\partial a} =
    (1-w)^{1-b}w^{1-a}
      \left(
        w^a\Gamma(a)^2 {}_3\tilde{F}_2(a,a,1-b;a+1,a+1;w)
        - B(a,b)I_w(a,b)\left(\log(w)-\psi(a) + \psi(a+b)\right)
      \right)/;w=I_z^{-1}(a,b)
   \f]
   \f[
   \frac{\partial }{\partial b} =
    (1-w)^{-b}w^{1-a}(w-1)
      \left(
        (1-w)^{b}\Gamma(b)^2 {}_3\tilde{F}_2(b,b,1-a;b+1,b+1;1-w)
        - B_{1-w}(b,a)\left(\log(1-w)-\psi(b) + \psi(a+b)\right)
      \right)/;w=I_z^{-1}(a,b)
   \f]
   \f[
   \frac{\partial }{\partial z} = (1-w)^{1-b}w^{1-a}B(a,b)/;w=I_z^{-1}(a,b)
   \f]
 *
 * @param a Shape parameter a >= 0; a and b can't both be 0
 * @param b Shape parameter b >= 0
 * @param p Random variate. 0 <= p <= 1
 * @throws if constraints are violated or if any argument is NaN
 * @return The inverse of the normalized incomplete beta function.
 */
template <typename T1, typename T2, typename T3,
          require_all_stan_scalar_t<T1, T2, T3>* = nullptr,
          require_any_var_t<T1, T2, T3>* = nullptr>
inline var inv_inc_beta(const T1& a, const T2& b, const T3& p) {
  double a_val = value_of(a);
  double b_val = value_of(b);
  double p_val = value_of(p);
  double w = inv_inc_beta(a_val, b_val, p_val);
  return make_callback_var(w, [a, b, p, a_val, b_val, w](auto& vi) {
    double log_w = log(w);
    double log1m_w = log1m(w);
    double one_m_a = 1 - a_val;
    double one_m_b = 1 - b_val;
    double one_m_w = 1 - w;
    double ap1 = a_val + 1;
    double bp1 = b_val + 1;
    double lbeta_ab = lbeta(a_val, b_val);
    double digamma_apb = digamma(a_val + b_val);

    if (!is_constant_all<T1>::value) {
      double da1 = exp(one_m_b * log1m_w + one_m_a * log_w);
      double da2
          = a_val * log_w + 2 * lgamma(a_val)
            + log(hypergeometric_3F2({a_val, a_val, one_m_b}, {ap1, ap1}, w))
            - 2 * lgamma(ap1);
      double da3 = inc_beta(a_val, b_val, w) * exp(lbeta_ab)
                   * (log_w - digamma(a_val) + digamma_apb);

      forward_as<var>(a).adj() += vi.adj() * da1 * (exp(da2) - da3);
    }

    if (!is_constant_all<T2>::value) {
      double db1 = (w - 1) * exp(-b_val * log1m_w + one_m_a * log_w);
      double db2 = 2 * lgamma(b_val)
                   + log(hypergeometric_3F2({b_val, b_val, one_m_a}, {bp1, bp1},
                                            one_m_w))
                   - 2 * lgamma(bp1) + b_val * log1m_w;

      double db3 = inc_beta(b_val, a_val, one_m_w) * exp(lbeta_ab)
                   * (log1m_w - digamma(b_val) + digamma_apb);

      forward_as<var>(b).adj() += vi.adj() * db1 * (exp(db2) - db3);
    }

    if (!is_constant_all<T3>::value) {
      forward_as<var>(p).adj()
          += vi.adj() * exp(one_m_b * log1m_w + one_m_a * log_w + lbeta_ab);
    }
  });
}

}  // namespace math
}  // namespace stan
#endif
