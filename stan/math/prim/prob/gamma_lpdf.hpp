#ifndef STAN_MATH_PRIM_PROB_GAMMA_LPDF_HPP
#define STAN_MATH_PRIM_PROB_GAMMA_LPDF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/as_column_vector_or_scalar.hpp>
#include <stan/math/prim/fun/as_array_or_scalar.hpp>
#include <stan/math/prim/fun/as_value_column_array_or_scalar.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/digamma.hpp>
#include <stan/math/prim/fun/grad_reg_inc_gamma.hpp>
#include <stan/math/prim/fun/lgamma.hpp>
#include <stan/math/prim/fun/log.hpp>
#include <stan/math/prim/fun/max_size.hpp>
#include <stan/math/prim/fun/scalar_seq_view.hpp>
#include <stan/math/prim/fun/size.hpp>
#include <stan/math/prim/fun/size_zero.hpp>
#include <stan/math/prim/fun/to_ref.hpp>
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
template <bool propto, typename T_y, typename T_shape, typename T_inv_scale,
          require_all_not_nonscalar_prim_or_rev_kernel_expression_t<
              T_y, T_shape, T_inv_scale>* = nullptr>
return_type_t<T_y, T_shape, T_inv_scale> gamma_lpdf(const T_y& y,
                                                    const T_shape& alpha,
                                                    const T_inv_scale& beta) {
  using T_partials_return = partials_return_t<T_y, T_shape, T_inv_scale>;
  using T_y_ref = ref_type_if_t<!is_constant<T_y>::value, T_y>;
  using T_alpha_ref = ref_type_if_t<!is_constant<T_shape>::value, T_shape>;
  using T_beta_ref
      = ref_type_if_t<!is_constant<T_inv_scale>::value, T_inv_scale>;
  static const char* function = "gamma_lpdf";
  check_consistent_sizes(function, "Random variable", y, "Shape parameter",
                         alpha, "Inverse scale parameter", beta);
  T_y_ref y_ref = y;
  T_alpha_ref alpha_ref = alpha;
  T_beta_ref beta_ref = beta;

  decltype(auto) y_val = to_ref(as_value_column_array_or_scalar(y_ref));
  decltype(auto) alpha_val = to_ref(as_value_column_array_or_scalar(alpha_ref));
  decltype(auto) beta_val = to_ref(as_value_column_array_or_scalar(beta_ref));

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
    logp = -sum(lgamma(alpha_val)) * N / math::size(alpha);
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
