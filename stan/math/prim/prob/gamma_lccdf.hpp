#ifndef STAN_MATH_PRIM_PROB_GAMMA_LCCDF_HPP
#define STAN_MATH_PRIM_PROB_GAMMA_LCCDF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/digamma.hpp>
#include <stan/math/prim/fun/exp.hpp>
#include <stan/math/prim/fun/gamma_q.hpp>
#include <stan/math/prim/fun/grad_reg_inc_gamma.hpp>
#include <stan/math/prim/fun/log.hpp>
#include <stan/math/prim/fun/max_size.hpp>
#include <stan/math/prim/fun/scalar_seq_view.hpp>
#include <stan/math/prim/fun/size.hpp>
#include <stan/math/prim/fun/size_zero.hpp>
#include <stan/math/prim/fun/tgamma.hpp>
#include <stan/math/prim/fun/value_of.hpp>
#include <stan/math/prim/functor/partials_propagator.hpp>
#include <cmath>

namespace stan {
namespace math {

template <typename T_y, typename T_shape, typename T_inv_scale>
inline return_type_t<T_y, T_shape, T_inv_scale> gamma_lccdf(
    const T_y& y, const T_shape& alpha, const T_inv_scale& beta) {
  using T_partials_return = partials_return_t<T_y, T_shape, T_inv_scale>;
  using std::exp;
  using std::log;
  using std::pow;
  using T_y_ref = ref_type_t<T_y>;
  using T_alpha_ref = ref_type_t<T_shape>;
  using T_beta_ref = ref_type_t<T_inv_scale>;
  static constexpr const char* function = "gamma_lccdf";
  check_consistent_sizes(function, "Random variable", y, "Shape parameter",
                         alpha, "Inverse scale parameter", beta);
  T_y_ref y_ref = y;
  T_alpha_ref alpha_ref = alpha;
  T_beta_ref beta_ref = beta;
  check_positive_finite(function, "Shape parameter", alpha_ref);
  check_positive_finite(function, "Inverse scale parameter", beta_ref);
  check_nonnegative(function, "Random variable", y_ref);

  if (size_zero(y, alpha, beta)) {
    return 0;
  }

  T_partials_return P(0.0);
  auto ops_partials = make_partials_propagator(y_ref, alpha_ref, beta_ref);

  scalar_seq_view<T_y_ref> y_vec(y_ref);
  scalar_seq_view<T_alpha_ref> alpha_vec(alpha_ref);
  scalar_seq_view<T_beta_ref> beta_vec(beta_ref);
  size_t N = max_size(y, alpha, beta);

  // Explicit return for extreme values
  // The gradients are technically ill-defined, but treated as zero
  for (size_t i = 0; i < stan::math::size(y); i++) {
    if (y_vec.val(i) == 0) {
      // LCCDF(0) = log(P(Y > 0)) = log(1) = 0
      return ops_partials.build(0.0);
    }
  }

  for (size_t n = 0; n < N; n++) {
    // Explicit results for extreme values
    // The gradients are technically ill-defined, but treated as zero
    if (y_vec.val(n) == INFTY) {
      // LCCDF(∞) = log(P(Y > ∞)) = log(0) = -∞
      return ops_partials.build(negative_infinity());
    }

    const T_partials_return y_dbl = y_vec.val(n);
    const T_partials_return alpha_dbl = alpha_vec.val(n);
    const T_partials_return beta_dbl = beta_vec.val(n);
    const T_partials_return beta_y_dbl = beta_dbl * y_dbl;

    // Qn = 1 - Pn
    const T_partials_return Qn = gamma_q(alpha_dbl, beta_y_dbl);
    const T_partials_return log_Qn = log(Qn);

    P += log_Qn;

    if constexpr (is_any_autodiff_v<T_y, T_inv_scale>) {
      const T_partials_return log_y_dbl = log(y_dbl);
      const T_partials_return log_beta_dbl = log(beta_dbl);
      const T_partials_return log_pdf
          = alpha_dbl * log_beta_dbl - lgamma(alpha_dbl)
            + (alpha_dbl - 1.0) * log_y_dbl - beta_y_dbl;
      const T_partials_return common_term = exp(log_pdf - log_Qn);

      if constexpr (is_autodiff_v<T_y>) {
        // d/dy log(1-F(y)) = -f(y)/(1-F(y))
        partials<0>(ops_partials)[n] -= common_term;
      }
      if constexpr (is_autodiff_v<T_inv_scale>) {
        // d/dbeta log(1-F(y)) = -y*f(y)/(beta*(1-F(y)))
        partials<2>(ops_partials)[n] -= y_dbl / beta_dbl * common_term;
      }
    }

    if constexpr (is_autodiff_v<T_shape>) {
      const T_partials_return digamma_val = digamma(alpha_dbl);
      const T_partials_return gamma_val = tgamma(alpha_dbl);
      // d/dalpha log(1-F(y)) = grad_upper_inc_gamma / (1-F(y))
      partials<1>(ops_partials)[n]
          += grad_reg_inc_gamma(alpha_dbl, beta_y_dbl, gamma_val, digamma_val)
             / Qn;
    }
  }
  return ops_partials.build(P);
}

}  // namespace math
}  // namespace stan

#endif
