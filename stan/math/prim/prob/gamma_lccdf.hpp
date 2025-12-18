#ifndef STAN_MATH_PRIM_PROB_GAMMA_LCCDF_HPP
#define STAN_MATH_PRIM_PROB_GAMMA_LCCDF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/digamma.hpp>
#include <stan/math/prim/fun/exp.hpp>
#include <stan/math/prim/fun/gamma_p.hpp>
#include <stan/math/prim/fun/grad_reg_lower_inc_gamma.hpp>
#include <stan/math/prim/fun/log.hpp>
#include <stan/math/prim/fun/log1m.hpp>
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
    return 0.0;
  }

  T_partials_return P(0.0);
  auto ops_partials = make_partials_propagator(y_ref, alpha_ref, beta_ref);

  scalar_seq_view<T_y_ref> y_vec(y_ref);
  scalar_seq_view<T_alpha_ref> alpha_vec(alpha_ref);
  scalar_seq_view<T_beta_ref> beta_vec(beta_ref);
  const size_t N = max_size(y, alpha, beta);

  for (size_t i = 0; i < stan::math::size(y); ++i) {
    if (y_vec.val(i) == 0.0) {
      return ops_partials.build(0.0);
    }
  }

  constexpr bool need_y_beta_deriv
      = is_any_autodiff_v<T_y, T_inv_scale>;
  constexpr bool alpha_is_scalar = is_constant_all<T_shape>::value;

  T_partials_return lgamma_alpha_const = 0.0;
  if constexpr (need_y_beta_deriv && alpha_is_scalar) {
    const T_partials_return alpha0 = value_of(alpha_vec.val(0));
    lgamma_alpha_const = lgamma(alpha0);
  }

  for (size_t n = 0; n < N; ++n) {
    const T_partials_return y_dbl = value_of(y_vec.val(n));

    if (y_dbl == INFTY) {
      return ops_partials.build(negative_infinity());
    }

    const T_partials_return alpha_dbl = value_of(alpha_vec.val(n));
    const T_partials_return beta_dbl = value_of(beta_vec.val(n));
    const T_partials_return beta_y = beta_dbl * y_dbl;

    // ---------- 1. VALUE: log CCDF via lower regularized gamma ----------
    // Pn = P(alpha, beta*y) = CDF
    // Qn = 1 - Pn = CCDF
    const T_partials_return Pn = gamma_p(alpha_dbl, beta_y);
    const T_partials_return log_Qn = log1m(Pn);   // = log(1 - Pn)
    const T_partials_return Qn = 1.0 - Pn;        // needed for gradients

    // If Qn underflows to 0 numerically, the log-CCDF is -infinity
    if (Qn <= 0.0) {
      return ops_partials.build(negative_infinity());
    }

    P += log_Qn;

    if constexpr (need_y_beta_deriv) {
      const T_partials_return log_y = log(y_dbl);
      const T_partials_return log_beta = log(beta_dbl);
      const T_partials_return lgamma_alpha
          = (alpha_is_scalar ? lgamma_alpha_const : lgamma(alpha_dbl));

      const T_partials_return log_pdf
          = alpha_dbl * log_beta - lgamma_alpha
            + (alpha_dbl - 1.0) * log_y - beta_y;

      // hazard = f(y) / Q(y), on the log scale as exp(log_pdf - log_Qn)
      const T_partials_return hazard = exp(log_pdf - log_Qn);

      if constexpr (is_autodiff_v<T_y>) {
        partials<0>(ops_partials)[n] += -hazard;
      }
      if constexpr (is_autodiff_v<T_inv_scale>) {
        partials<2>(ops_partials)[n] += -(y_dbl / beta_dbl) * hazard;
      }
    }

    if constexpr (is_autodiff_v<T_shape>) {
      // For the shape derivative, we stay entirely on the P side:
      //
      //   Q(alpha, z) = 1 - P(alpha, z)
      //   d/da Q = - d/da P
      //
      // so
      //
      //   d/da log Q = (1 / Q) * dQ/da
      //              = - (1 / Q) * dP/da.
      //
      // grad_reg_lower_inc_gamma(alpha, z) = d/da P(alpha, z)
      const T_partials_return dP_da
          = grad_reg_lower_inc_gamma(alpha_dbl, beta_y);
      partials<1>(ops_partials)[n] -= dP_da / Qn;
    }
  }
  return ops_partials.build(P);
}

}  // namespace math
}  // namespace stan

#endif