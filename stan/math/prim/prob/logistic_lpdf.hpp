#ifndef STAN_MATH_PRIM_PROB_LOGISTIC_LPDF_HPP
#define STAN_MATH_PRIM_PROB_LOGISTIC_LPDF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/as_column_vector_or_scalar.hpp>
#include <stan/math/prim/fun/as_array_or_scalar.hpp>
#include <stan/math/prim/fun/as_value_column_array_or_scalar.hpp>
#include <stan/math/prim/fun/exp.hpp>
#include <stan/math/prim/fun/log.hpp>
#include <stan/math/prim/fun/log1p.hpp>
#include <stan/math/prim/fun/log1p_exp.hpp>
#include <stan/math/prim/fun/max_size.hpp>
#include <stan/math/prim/fun/size.hpp>
#include <stan/math/prim/fun/size_zero.hpp>
#include <stan/math/prim/fun/to_ref.hpp>
#include <stan/math/prim/fun/value_of.hpp>
#include <stan/math/prim/functor/partials_propagator.hpp>
#include <cmath>

namespace stan {
namespace math {

// Logistic(y|mu, sigma)    [sigma > 0]
template <bool propto, typename T_y, typename T_loc, typename T_scale,
          require_all_not_nonscalar_prim_or_rev_kernel_expression_t<
              T_y, T_loc, T_scale>* = nullptr>
return_type_t<T_y, T_loc, T_scale> logistic_lpdf(const T_y& y, const T_loc& mu,
                                                 const T_scale& sigma) {
  using T_partials_return = partials_return_t<T_y, T_loc, T_scale>;
  using T_y_ref = ref_type_if_not_constant_t<T_y>;
  using T_mu_ref = ref_type_if_not_constant_t<T_loc>;
  using T_sigma_ref = ref_type_if_not_constant_t<T_scale>;
  static constexpr const char* function = "logistic_lpdf";
  check_consistent_sizes(function, "Random variable", y, "Location parameter",
                         mu, "Scale parameter", sigma);

  T_y_ref y_ref = y;
  T_mu_ref mu_ref = mu;
  T_sigma_ref sigma_ref = sigma;

  decltype(auto) y_val = to_ref(as_value_column_array_or_scalar(y_ref));
  decltype(auto) mu_val = to_ref(as_value_column_array_or_scalar(mu_ref));
  decltype(auto) sigma_val = to_ref(as_value_column_array_or_scalar(sigma_ref));

  check_finite(function, "Random variable", y_val);
  check_finite(function, "Location parameter", mu_val);
  check_positive_finite(function, "Scale parameter", sigma_val);

  if (size_zero(y, mu, sigma)) {
    return 0.0;
  }
  if (!include_summand<propto, T_y, T_loc, T_scale>::value) {
    return 0.0;
  }

  auto ops_partials = make_partials_propagator(y_ref, mu_ref, sigma_ref);

  const auto& inv_sigma
      = to_ref_if<!is_constant_all<T_y, T_loc, T_scale>::value>(inv(sigma_val));
  const auto& y_minus_mu
      = to_ref_if<!is_constant_all<T_scale>::value>(y_val - mu_val);
  const auto& y_minus_mu_div_sigma = to_ref(y_minus_mu * inv_sigma);

  size_t N = max_size(y, mu, sigma);
  T_partials_return logp = -sum(y_minus_mu_div_sigma)
                           - 2.0 * sum(log1p_exp(-y_minus_mu_div_sigma));
  if (include_summand<propto, T_scale>::value) {
    logp -= sum(log(sigma_val)) * N / math::size(sigma);
  }

  if (!is_constant_all<T_y, T_scale>::value) {
    const auto& exp_y_minus_mu_div_sigma = exp(y_minus_mu_div_sigma);
    const auto& y_deriv = to_ref_if<(!is_constant_all<T_scale>::value
                                     && !is_constant_all<T_y>::value)>(
        (2 / (1 + exp_y_minus_mu_div_sigma) - 1) * inv_sigma);
    if (!is_constant_all<T_y>::value) {
      partials<0>(ops_partials) = y_deriv;
    }
    if (!is_constant_all<T_scale>::value) {
      partials<2>(ops_partials) = (-y_deriv * y_minus_mu - 1) * inv_sigma;
    }
  }
  if (!is_constant_all<T_loc>::value) {
    const auto& exp_mu_div_sigma = to_ref(exp(mu_val * inv_sigma));
    edge<1>(ops_partials).partials_
        = (1
           - 2 * exp_mu_div_sigma / (exp_mu_div_sigma + exp(y_val * inv_sigma)))
          * inv_sigma;
  }
  return ops_partials.build(logp);
}

template <typename T_y, typename T_loc, typename T_scale>
inline return_type_t<T_y, T_loc, T_scale> logistic_lpdf(const T_y& y,
                                                        const T_loc& mu,
                                                        const T_scale& sigma) {
  return logistic_lpdf<false>(y, mu, sigma);
}

}  // namespace math
}  // namespace stan
#endif
