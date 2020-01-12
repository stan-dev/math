#ifndef STAN_MATH_PRIM_PROB_GAMMA_LPDF_HPP
#define STAN_MATH_PRIM_PROB_GAMMA_LPDF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/size_zero.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/value_of.hpp>
#include <stan/math/prim/fun/digamma.hpp>
#include <stan/math/prim/fun/lgamma.hpp>
#include <stan/math/prim/fun/grad_reg_inc_gamma.hpp>
#include <cmath>

namespace stan {
namespace math {

/** \ingroup prob_dists
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
return_type_t<T_y, T_shape, T_inv_scale> gamma_lpdf(const T_y& y,
                                                    const T_shape& alpha,
                                                    const T_inv_scale& beta) {
  static const char* function = "gamma_lpdf";
  using T_partials_return = partials_return_t<T_y, T_shape, T_inv_scale>;

  if (size_zero(y, alpha, beta)) {
    return 0.0;
  }

  T_partials_return logp(0.0);

  check_not_nan(function, "Random variable", y);
  check_positive_finite(function, "Shape parameter", alpha);
  check_positive_finite(function, "Inverse scale parameter", beta);
  check_consistent_sizes(function, "Random variable", y, "Shape parameter",
                         alpha, "Inverse scale parameter", beta);

  if (!include_summand<propto, T_y, T_shape, T_inv_scale>::value) {
    return 0.0;
  }

  scalar_seq_view<T_y> y_vec(y);
  scalar_seq_view<T_shape> alpha_vec(alpha);
  scalar_seq_view<T_inv_scale> beta_vec(beta);

  for (size_t n = 0; n < size(y); n++) {
    const T_partials_return y_dbl = value_of(y_vec[n]);
    if (y_dbl < 0) {
      return LOG_ZERO;
    }
  }

  size_t N = max_size(y, alpha, beta);
  operands_and_partials<T_y, T_shape, T_inv_scale> ops_partials(y, alpha, beta);

  using std::log;

  VectorBuilder<include_summand<propto, T_y, T_shape>::value, T_partials_return,
                T_y>
      log_y(size(y));
  if (include_summand<propto, T_y, T_shape>::value) {
    for (size_t n = 0; n < size(y); n++) {
      if (value_of(y_vec[n]) > 0) {
        log_y[n] = log(value_of(y_vec[n]));
      }
    }
  }

  VectorBuilder<include_summand<propto, T_shape>::value, T_partials_return,
                T_shape>
      lgamma_alpha(size(alpha));
  VectorBuilder<!is_constant_all<T_shape>::value, T_partials_return, T_shape>
      digamma_alpha(size(alpha));
  for (size_t n = 0; n < size(alpha); n++) {
    if (include_summand<propto, T_shape>::value) {
      lgamma_alpha[n] = lgamma(value_of(alpha_vec[n]));
    }
    if (!is_constant_all<T_shape>::value) {
      digamma_alpha[n] = digamma(value_of(alpha_vec[n]));
    }
  }

  VectorBuilder<include_summand<propto, T_shape, T_inv_scale>::value,
                T_partials_return, T_inv_scale>
      log_beta(size(beta));
  if (include_summand<propto, T_shape, T_inv_scale>::value) {
    for (size_t n = 0; n < size(beta); n++) {
      log_beta[n] = log(value_of(beta_vec[n]));
    }
  }

  for (size_t n = 0; n < N; n++) {
    const T_partials_return y_dbl = value_of(y_vec[n]);
    const T_partials_return alpha_dbl = value_of(alpha_vec[n]);
    const T_partials_return beta_dbl = value_of(beta_vec[n]);

    if (include_summand<propto, T_shape>::value) {
      logp -= lgamma_alpha[n];
    }
    if (include_summand<propto, T_shape, T_inv_scale>::value) {
      logp += alpha_dbl * log_beta[n];
    }
    if (include_summand<propto, T_y, T_shape>::value) {
      logp += (alpha_dbl - 1.0) * log_y[n];
    }
    if (include_summand<propto, T_y, T_inv_scale>::value) {
      logp -= beta_dbl * y_dbl;
    }

    if (!is_constant_all<T_y>::value) {
      ops_partials.edge1_.partials_[n] += (alpha_dbl - 1) / y_dbl - beta_dbl;
    }
    if (!is_constant_all<T_shape>::value) {
      ops_partials.edge2_.partials_[n]
          += -digamma_alpha[n] + log_beta[n] + log_y[n];
    }
    if (!is_constant_all<T_inv_scale>::value) {
      ops_partials.edge3_.partials_[n] += alpha_dbl / beta_dbl - y_dbl;
    }
  }
  return ops_partials.build(logp);
}

template <typename T_y, typename T_shape, typename T_inv_scale>
inline return_type_t<T_y, T_shape, T_inv_scale> gamma_lpdf(
    const T_y& y, const T_shape& alpha, const T_inv_scale& beta) {
  return gamma_lpdf<false>(y, alpha, beta);
}

}  // namespace math
}  // namespace stan
#endif
