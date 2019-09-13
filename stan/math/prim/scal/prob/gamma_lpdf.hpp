#ifndef STAN_MATH_PRIM_SCAL_PROB_GAMMA_LPDF_HPP
#define STAN_MATH_PRIM_SCAL_PROB_GAMMA_LPDF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/scal/err/check_consistent_sizes.hpp>
#include <stan/math/prim/scal/err/check_not_nan.hpp>
#include <stan/math/prim/scal/err/check_positive_finite.hpp>
#include <stan/math/prim/scal/fun/size_zero.hpp>
#include <stan/math/prim/scal/fun/constants.hpp>
#include <stan/math/prim/scal/fun/value_of.hpp>
#include <stan/math/prim/scal/fun/digamma.hpp>
#include <stan/math/prim/scal/fun/lgamma.hpp>
#include <stan/math/prim/scal/fun/grad_reg_inc_gamma.hpp>
#include <cmath>

namespace stan {
namespace math {

/**
 * The log of a gamma density for y with the specified
 * shape and inverse scale parameters.
 * Shape and inverse scale parameters must be greater than 0.
 * y must be greater than or equal to 0.
 *
 \f{eqnarray*}{
 y &\sim& \mbox{\sf{Gamma}}(\alpha, \beta) \\
 \log (p (y \, |\, \alpha, \beta) ) &=& \log \left(
 \frac{\beta^\alpha}{\Gamma(\alpha)} y^{\alpha - 1} \exp^{- \beta y} \right) \\
 &=& \alpha \log(\beta) - \log(\Gamma(\alpha)) + (\alpha - 1) \log(y) - \beta
 y\\ & & \mathrm{where} \; y > 0 \f}
 * @param y A scalar variable.
 * @param alpha Shape parameter.
 * @param beta Inverse scale parameter.
 * @throw std::domain_error if alpha is not greater than 0.
 * @throw std::domain_error if beta is not greater than 0.
 * @throw std::domain_error if y is not greater than or equal to 0.
 * @tparam T_y Type of scalar.
 * @tparam T_shape Type of shape.
 * @tparam T_inv_scale Type of inverse scale.
 */
template <bool propto, typename T_y, typename T_shape, typename T_inv_scale>
inline auto gamma_lpdf(const T_y& y, const T_shape& alpha,
                       const T_inv_scale& beta) {
  using T_partials = partials_return_t<T_y, T_shape, T_inv_scale>;
  T_partials logp(0.0);

  using std::log;

  static const char* function = "gamma_lpdf";
  check_not_nan(function, "Random variable", y);
  check_positive_finite(function, "Shape parameter", alpha);
  check_positive_finite(function, "Inverse scale parameter", beta);
  check_consistent_sizes(function, "Random variable", y, "Shape parameter",
                         alpha, "Inverse scale parameter", beta);

  const scalar_seq_view<T_y> y_vec(y);
  const scalar_seq_view<T_shape> alpha_vec(alpha);
  const scalar_seq_view<T_inv_scale> beta_vec(beta);
  const size_t N = max_size(y, alpha, beta);
  const size_t size_y = length(y);
  operands_and_partials<T_y, T_shape, T_inv_scale> ops_partials(y, alpha, beta);
  if (!include_summand<propto, T_y, T_shape, T_inv_scale>::value) {
    return ops_partials.build(logp);
  } else if (size_zero(y, alpha, beta)) {
    return ops_partials.build(logp);
  }
  T_partials log_beta = 0;
  for (size_t n = 0; n < N; n++) {
    const T_partials y_dbl = value_of(y_vec[n]);
    const T_partials alpha_dbl = value_of(alpha_vec[n]);
    const T_partials beta_dbl = value_of(beta_vec[n]);
    T_partials log_y = 0;
    if (y_dbl < 0) {
      return ops_partials.build(T_partials(LOG_ZERO));
    }
    if (include_summand<propto, T_shape>::value) {
      const T_partials lgamma_alpha = lgamma(alpha_dbl);
      logp -= lgamma_alpha;
    }
    if (include_summand<propto, T_shape, T_inv_scale>::value) {
      log_beta = log(beta_dbl);
      logp += alpha_dbl * log_beta;
    }
    if (include_summand<propto, T_y, T_shape>::value) {
      if (y_dbl > 0) {
        log_y = log(y_dbl);
      }
      logp += (alpha_dbl - 1.0) * log_y;
    }
    if (include_summand<propto, T_y, T_inv_scale>::value) {
      logp -= beta_dbl * y_dbl;
    }

    if (!is_constant_all<T_y>::value) {
      ops_partials.edge1_.partials_[n] += (alpha_dbl - 1) / y_dbl - beta_dbl;
    }
    if (!is_constant_all<T_shape>::value) {
      const T_partials digamma_alpha = digamma(alpha_dbl);
      ops_partials.edge2_.partials_[n] += -digamma_alpha + log_beta + log_y;
    }
    if (!is_constant_all<T_inv_scale>::value) {
      ops_partials.edge3_.partials_[n] += alpha_dbl / beta_dbl - y_dbl;
    }
  }
  return ops_partials.build(logp);
}

template <typename T_y, typename T_shape, typename T_inv_scale>
inline auto gamma_lpdf(const T_y& y, const T_shape& alpha,
                       const T_inv_scale& beta) {
  return gamma_lpdf<false>(y, alpha, beta);
}

}  // namespace math
}  // namespace stan
#endif
