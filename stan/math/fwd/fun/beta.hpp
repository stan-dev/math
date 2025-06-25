#ifndef STAN_MATH_FWD_FUN_BETA_HPP
#define STAN_MATH_FWD_FUN_BETA_HPP

#include <stan/math/fwd/meta.hpp>
#include <stan/math/fwd/core.hpp>
#include <stan/math/fwd/fun/digamma.hpp>
#include <stan/math/prim/fun/beta.hpp>

namespace stan {
namespace math {

/**
 * Return fvar with the beta function applied to the specified
 * arguments and its gradient.
 *
 * The beta function is defined for \f$a > 0\f$ and \f$b > 0\f$ by
 *
 * \f$\mbox{B}(a, b) = \frac{\Gamma(a) \Gamma(b)}{\Gamma(a+b)}\f$.
 *
   \f[
   \mbox{beta}(\alpha, \beta) =
   \begin{cases}
     \int_0^1 u^{\alpha - 1} (1 - u)^{\beta - 1} \, du & \mbox{if } \alpha,
 \beta>0 \\[6pt] \textrm{NaN} & \mbox{if } \alpha = \textrm{NaN or } \beta =
 \textrm{NaN} \end{cases} \f]

   \f[
   \frac{\partial\, \mbox{beta}(\alpha, \beta)}{\partial \alpha} =
   \begin{cases}
     \left(\psi(\alpha)-\psi(\alpha+\beta)\right)*\mbox{beta}(\alpha, \beta)
      & \mbox{if } \alpha, \beta>0 \\[6pt]
     \textrm{NaN} & \mbox{if } \alpha = \textrm{NaN or } \beta = \textrm{NaN}
   \end{cases}
   \f]

   \f[
   \frac{\partial\, \mbox{beta}(\alpha, \beta)}{\partial \beta} =
   \begin{cases}
     \left(\psi(\beta)-\psi(\alpha+\beta)\right)*\mbox{beta}(\alpha, \beta)
      & \mbox{if } \alpha, \beta>0 \\[6pt]
     \textrm{NaN} & \mbox{if } \alpha = \textrm{NaN or } \beta = \textrm{NaN}
   \end{cases}
   \f]
 *
 * @tparam T inner type of the fvar
 * @param x1 First value
 * @param x2 Second value
 * @return Fvar with result beta function of arguments and gradients.
 */
template <typename Ta, typename Tb,
          typename FvarInnerT = partials_return_t<Ta, Tb>,
          require_return_type_t<is_fvar, Ta, Tb>* = nullptr,
          require_all_stan_scalar_t<Ta, Tb>* = nullptr>
inline fvar<FvarInnerT> beta(const Ta& a, const Tb& b) {
  const auto& a_val = value_of(a);
  const auto& b_val = value_of(b);
  const FvarInnerT beta_val = beta(a_val, b_val);
  const FvarInnerT digamma_ab = digamma(a_val + b_val);
  FvarInnerT beta_d(0);
  if constexpr (!is_constant<Ta>::value) {
    beta_d += (digamma(a_val) - digamma_ab) * beta_val * a.d_;
  }
  if constexpr (!is_constant<Tb>::value) {
    beta_d += (digamma(b_val) - digamma_ab) * beta_val * b.d_;
  }
  return fvar<FvarInnerT>(beta_val, beta_d);
}

}  // namespace math
}  // namespace stan
#endif
