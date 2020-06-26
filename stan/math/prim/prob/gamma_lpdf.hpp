#ifndef STAN_MATH_PRIM_PROB_GAMMA_LPDF_HPP
#define STAN_MATH_PRIM_PROB_GAMMA_LPDF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/digamma.hpp>
#include <stan/math/prim/fun/grad_reg_inc_gamma.hpp>
#include <stan/math/prim/fun/lgamma.hpp>
#include <stan/math/prim/fun/log.hpp>
#include <stan/math/prim/fun/max_size.hpp>
#include <stan/math/prim/fun/size.hpp>
#include <stan/math/prim/fun/size_zero.hpp>
#include <stan/math/prim/fun/value_of.hpp>
#include <stan/math/prim/functor/operands_and_partials.hpp>
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
 *
 * @tparam T_y type of scalar
 * @tparam T_shape type of shape
 * @tparam T_inv_scale type of inverse scale
 * @param y A scalar variable.
 * @param alpha Shape parameter.
 * @param beta Inverse scale parameter.
 * @throw std::domain_error if alpha is not greater than 0.
 * @throw std::domain_error if beta is not greater than 0.
 * @throw std::domain_error if y is not greater than or equal to 0.
 */
template <bool propto, typename T_y, typename T_shape, typename T_inv_scale>
return_type_t<T_y, T_shape, T_inv_scale> gamma_lpdf(const T_y& y,
                                                    const T_shape& alpha,
                                                    const T_inv_scale& beta) {
  using T_partials_return = partials_return_t<T_y, T_shape, T_inv_scale>;
  using T_y_ref = ref_type_if_t<!is_constant<T_y>::value, T_y>;
  using T_alpha_ref
      = ref_type_if_t<!is_constant<T_shape>::value, T_shape>;
  using T_beta_ref
      = ref_type_if_t<!is_constant<T_inv_scale>::value, T_inv_scale>;
  static const char* function = "gamma_lpdf";
  check_consistent_sizes(function, "Random variable", y, "Shape parameter",
                         alpha, "Inverse scale parameter", beta);
  T_y_ref y_ref = y;
  T_alpha_ref alpha_ref = alpha;
  T_beta_ref beta_ref = beta;

  const auto& y_col = as_column_vector_or_scalar(y_ref);
  const auto& alpha_col = as_column_vector_or_scalar(alpha_ref);
  const auto& beta_col = as_column_vector_or_scalar(beta_ref);

  const auto& y_arr = as_array_or_scalar(y_col);
  const auto& alpha_arr = as_array_or_scalar(alpha_col);
  const auto& beta_arr = as_array_or_scalar(beta_col);

  ref_type_t<decltype(value_of(y_arr))> y_val = value_of(y_arr);
  ref_type_t<decltype(value_of(alpha_arr))> alpha_val = value_of(alpha_arr);
  ref_type_t<decltype(value_of(beta_arr))> beta_val = value_of(beta_arr);

  check_not_nan(function, "Random variable", y_val);
  check_positive_finite(function, "Shape parameter", alpha_val);
  check_positive_finite(function, "Inverse scale parameter", beta_val);

  if (size_zero(y, alpha, beta)) {
    return 0.0;
  }
  if (!include_summand<propto, T_y, T_shape, T_inv_scale>::value) {
    return 0.0;
  }

  operands_and_partials<T_y_ref, T_alpha_ref, T_beta_ref> ops_partials(
      y_ref, alpha_ref, beta_ref);

  scalar_seq_view<decltype(y_val)> y_vec(y_val);
  for (size_t n = 0; n < stan::math::size(y); n++) {
    if (y_vec[n] < 0) {
      return LOG_ZERO;
    }
  }

  size_t N = max_size(y, alpha, beta);
  T_partials_return logp(0.0);
  if (include_summand<propto, T_shape>::value) {
    logp = -sum(lgamma(alpha_val)) * N / size(alpha);
  }
  const auto& log_y = to_ref_if<is_constant_all<T_shape>::value>(log(y_val));
  if (include_summand<propto, T_shape, T_inv_scale>::value) {
    const auto& log_beta
        = to_ref_if<!is_constant_all<T_shape>::value>(log(beta_val));
    logp += sum(alpha_val * log_beta) * N / max_size(alpha, beta);
    if (!is_constant_all<T_shape>::value) {
      ops_partials.edge2_.partials_ = log_beta + log_y - digamma(alpha_val);
    }
  }
  if (include_summand<propto, T_y, T_shape>::value) {
    logp += sum((alpha_val - 1.0) * log_y) * N / max_size(alpha, y);
  }
  if (include_summand<propto, T_y, T_inv_scale>::value) {
    logp -= sum(beta_val * y_val) * N / max_size(beta, y);
  }

  if (!is_constant_all<T_y>::value) {
    ops_partials.edge1_.partials_ = (alpha_val - 1) / y_val - beta_val;
  }
  if (!is_constant_all<T_inv_scale>::value) {
    ops_partials.edge3_.partials_ = alpha_val / beta_val - y_val;
  }

//  scalar_seq_view<T_y> y_vec(y);
//  scalar_seq_view<T_shape> alpha_vec(alpha);
//  scalar_seq_view<T_inv_scale> beta_vec(beta);
//  size_t N = max_size(y, alpha, beta);

//  for (size_t n = 0; n < stan::math::size(y); n++) {
//    const T_partials_return y_dbl = value_of(y_vec[n]);
//    if (y_dbl < 0) {
//      return LOG_ZERO;
//    }
//  }

//  VectorBuilder<include_summand<propto, T_y, T_shape>::value, T_partials_return,
//                T_y>
//      log_y(size(y));
//  if (include_summand<propto, T_y, T_shape>::value) {
//    for (size_t n = 0; n < stan::math::size(y); n++) {
//      if (value_of(y_vec[n]) > 0) {
//        log_y[n] = log(value_of(y_vec[n]));
//      }
//    }
//  }

//  VectorBuilder<include_summand<propto, T_shape>::value, T_partials_return,
//                T_shape>
//      lgamma_alpha(size(alpha));
//  VectorBuilder<!is_constant_all<T_shape>::value, T_partials_return, T_shape>
//      digamma_alpha(size(alpha));
//  for (size_t n = 0; n < stan::math::size(alpha); n++) {
//    if (include_summand<propto, T_shape>::value) {
//      lgamma_alpha[n] = lgamma(value_of(alpha_vec[n]));
//    }
//    if (!is_constant_all<T_shape>::value) {
//      digamma_alpha[n] = digamma(value_of(alpha_vec[n]));
//    }
//  }

//  VectorBuilder<include_summand<propto, T_shape, T_inv_scale>::value,
//                T_partials_return, T_inv_scale>
//      log_beta(size(beta));
//  if (include_summand<propto, T_shape, T_inv_scale>::value) {
//    for (size_t n = 0; n < stan::math::size(beta); n++) {
//      log_beta[n] = log(value_of(beta_vec[n]));
//    }
//  }

//  for (size_t n = 0; n < N; n++) {
//    const T_partials_return y_dbl = value_of(y_vec[n]);
//    const T_partials_return alpha_dbl = value_of(alpha_vec[n]);
//    const T_partials_return beta_dbl = value_of(beta_vec[n]);

//    if (include_summand<propto, T_shape>::value) {
//      logp -= lgamma_alpha[n];
//    }
//    if (include_summand<propto, T_shape, T_inv_scale>::value) {
//      logp += alpha_dbl * log_beta[n];
//    }
//    if (include_summand<propto, T_y, T_shape>::value) {
//      logp += (alpha_dbl - 1.0) * log_y[n];
//    }
//    if (include_summand<propto, T_y, T_inv_scale>::value) {
//      logp -= beta_dbl * y_dbl;
//    }

//    if (!is_constant_all<T_y>::value) {
//      ops_partials.edge1_.partials_[n] += (alpha_dbl - 1) / y_dbl - beta_dbl;
//    }
//    if (!is_constant_all<T_shape>::value) {
//      ops_partials.edge2_.partials_[n]
//          += -digamma_alpha[n] + log_beta[n] + log_y[n];
//    }
//    if (!is_constant_all<T_inv_scale>::value) {
//      ops_partials.edge3_.partials_[n] += alpha_dbl / beta_dbl - y_dbl;
//    }
//  }
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
