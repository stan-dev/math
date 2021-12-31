#ifndef STAN_MATH_FWD_FUN_INC_BETA_INV_HPP
#define STAN_MATH_FWD_FUN_INC_BETA_INV_HPP

#include <stan/math/fwd/meta.hpp>
#include <stan/math/fwd/core.hpp>
#include <stan/math/prim/fun/inc_beta_inv.hpp>
#include <stan/math/prim/fun/inc_beta.hpp>
#include <stan/math/prim/fun/exp.hpp>
#include <stan/math/prim/fun/log.hpp>
#include <stan/math/prim/fun/log_diff_exp.hpp>
#include <stan/math/prim/fun/lbeta.hpp>
#include <stan/math/prim/fun/lgamma.hpp>
#include <stan/math/prim/fun/digamma.hpp>
#include <stan/math/prim/fun/F32.hpp>

namespace stan {
namespace math {

/**
 * The fused multiply-add operation (C99).
 *
 * This double-based operation delegates to <code>fma</code>.
 *
 * The function is defined by
 *
 * <code>fma(a, b, c) = (a * b) + c</code>.
 *
 *
   \f[
   \mbox{fma}(x, y, z) =
   \begin{cases}
     x\cdot y+z & \mbox{if } -\infty\leq x, y, z \leq \infty \\[6pt]
     \textrm{NaN} & \mbox{if } x = \textrm{NaN}
   \end{cases}
   \f]

   \f[
   \frac{\partial\, \mbox{fma}(x, y, z)}{\partial x} =
   \begin{cases}
     y & \mbox{if } -\infty\leq x, y, z \leq \infty \\[6pt]
     \textrm{NaN} & \mbox{if } x = \textrm{NaN}
   \end{cases}
   \f]

   \f[
   \frac{\partial\, \mbox{fma}(x, y, z)}{\partial y} =
   \begin{cases}
     x & \mbox{if } -\infty\leq x, y, z \leq \infty \\[6pt]
     \textrm{NaN} & \mbox{if } x = \textrm{NaN}
   \end{cases}
   \f]

   \f[
   \frac{\partial\, \mbox{fma}(x, y, z)}{\partial z} =
   \begin{cases}
     1 & \mbox{if } -\infty\leq x, y, z \leq \infty \\[6pt]
     \textrm{NaN} & \mbox{if } x = \textrm{NaN}
   \end{cases}
   \f]
 *
 * @param x1 First value.
 * @param x2 Second value.
 * @param x3 Third value.
 * @return Product of the first two values plus the third.
 */
template <typename T1, typename T2, typename T3,
          require_all_stan_scalar_t<T1, T2, T3>* = nullptr,
          require_any_fvar_t<T1, T2, T3>* = nullptr>
inline fvar<partials_return_t<T1, T2, T3>> inc_beta_inv(const T1& a,
                                           const T2& b,
                                           const T3& p) {
  using T_return = partials_return_t<T1, T2, T3>;
  auto a_val = value_of(a);
  auto b_val = value_of(b);
  auto p_val = value_of(p);
  T_return w = inc_beta_inv(a_val, b_val, p_val);
  T_return log_w = log(w);
  T_return log1m_w = log1m(w);
  auto one_m_a = 1 - a_val;
  auto one_m_b = 1 - b_val;
  T_return one_m_w = 1 - w;
  auto ap1 = a_val + 1;
  auto bp1 = b_val + 1;
  auto lbeta_ab = lbeta(a_val, b_val);
  auto digamma_apb = digamma(a_val + b_val);

  T_return inv_d_(0);

  if(is_fvar<T1>::value) {
    auto da1 = exp(one_m_b * log1m_w + one_m_a * log_w);
    auto da2 = exp(a_val * log_w + 2 * lgamma(a_val)
                  + log(F32(a_val, a_val, one_m_b, ap1, ap1, w))
                  - 2 * lgamma(ap1));
    auto da3 =  inc_beta(a_val, b_val, w) * exp(lbeta_ab)
                  * (log_w - digamma(a_val) + digamma_apb);
    inv_d_ += forward_as<fvar<T_return>>(a).d_ * da1 * (da2 - da3);
  }

  if(is_fvar<T2>::value) {
    auto db1 = (w - 1) * exp(-b_val * log1m_w + one_m_a * log_w);
    auto db2 = 2 * lgamma(b_val) 
                  + log(F32(b_val, b_val, one_m_a, bp1, bp1, one_m_w))
                  - 2 * lgamma(bp1)
                  + b_val * log1m_w;

    auto db3 = inc_beta(b_val, a_val, one_m_w) * exp(lbeta_ab)
                  * (log1m_w - digamma(b_val) + digamma_apb);


    inv_d_ += forward_as<fvar<T_return>>(b).d_ * db1 * (exp(db2) - db3);
  }

  if(is_fvar<T3>::value) {
    inv_d_ += forward_as<fvar<T_return>>(p).d_ *
      exp(one_m_b * log1m_w + one_m_a * log_w + lbeta_ab);
  }

  return fvar<T_return>(w, inv_d_);
}


}  // namespace math
}  // namespace stan
#endif
