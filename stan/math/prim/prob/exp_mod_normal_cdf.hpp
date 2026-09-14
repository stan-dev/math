#ifndef STAN_MATH_PRIM_PROB_EXP_MOD_NORMAL_CDF_HPP
#define STAN_MATH_PRIM_PROB_EXP_MOD_NORMAL_CDF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/as_column_vector_or_scalar.hpp>
#include <stan/math/prim/fun/as_array_or_scalar.hpp>
#include <stan/math/prim/fun/as_value_column_array_or_scalar.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/erf.hpp>
#include <stan/math/prim/fun/exp.hpp>
#include <stan/math/prim/fun/inv.hpp>
#include <stan/math/prim/fun/is_inf.hpp>
#include <stan/math/prim/fun/log.hpp>
#include <stan/math/prim/fun/max_size.hpp>
#include <stan/math/prim/fun/scalar_seq_view.hpp>
#include <stan/math/prim/fun/size.hpp>
#include <stan/math/prim/fun/size_zero.hpp>
#include <stan/math/prim/fun/square.hpp>
#include <stan/math/prim/fun/to_ref.hpp>
#include <stan/math/prim/fun/value_of.hpp>
#include <stan/math/prim/fun/value_of_rec.hpp>
#include <stan/math/prim/functor/partials_propagator.hpp>
#include <stan/math/prim/prob/exp_mod_normal_utils.hpp>
#include <cmath>

namespace stan {
namespace math {

template <typename T_y, typename T_loc, typename T_scale, typename T_inv_scale,
          require_all_not_nonscalar_prim_or_rev_kernel_expression_t<
              T_y, T_loc, T_scale, T_inv_scale>* = nullptr>
inline return_type_t<T_y, T_loc, T_scale, T_inv_scale> exp_mod_normal_cdf(
    const T_y& y, const T_loc& mu, const T_scale& sigma,
    const T_inv_scale& lambda) {
  using T_partials_return = partials_return_t<T_y, T_loc, T_scale, T_inv_scale>;
  using std::log;
  using T_y_ref = ref_type_if_not_constant_t<T_y>;
  using T_mu_ref = ref_type_if_not_constant_t<T_loc>;
  using T_sigma_ref = ref_type_if_not_constant_t<T_scale>;
  using T_lambda_ref = ref_type_if_not_constant_t<T_inv_scale>;
  static constexpr const char* function = "exp_mod_normal_cdf";
  check_consistent_sizes(function, "Random variable", y, "Location parameter",
                         mu, "Scale parameter", sigma, "Inv_scale parameter",
                         lambda);
  T_y_ref y_ref = y;
  T_mu_ref mu_ref = mu;
  T_sigma_ref sigma_ref = sigma;
  T_lambda_ref lambda_ref = lambda;

  decltype(auto) y_val = to_ref(as_value_column_array_or_scalar(y_ref));
  decltype(auto) mu_val = to_ref(as_value_column_array_or_scalar(mu_ref));
  decltype(auto) sigma_val = to_ref(as_value_column_array_or_scalar(sigma_ref));
  decltype(auto) lambda_val
      = to_ref(as_value_column_array_or_scalar(lambda_ref));

  check_not_nan(function, "Random variable", y_val);
  check_finite(function, "Location parameter", mu_val);
  check_positive_finite(function, "Scale parameter", sigma_val);
  check_positive_finite(function, "Inv_scale parameter", lambda_val);

  if (size_zero(y, mu, sigma, lambda)) {
    return 1.0;
  }

  auto ops_partials
      = make_partials_propagator(y_ref, mu_ref, sigma_ref, lambda_ref);

  if constexpr (is_vector<T_y>::value) {
    if ((y_val == NEGATIVE_INFTY).any()) {
      return ops_partials.build(0.0);
    }
  } else {
    if (y_val == NEGATIVE_INFTY) {
      return ops_partials.build(0.0);
    }
  }

  scalar_seq_view<decltype(y_val)> y_vec(y_val);
  scalar_seq_view<decltype(mu_val)> mu_vec(mu_val);
  scalar_seq_view<decltype(sigma_val)> sigma_vec(sigma_val);
  scalar_seq_view<decltype(lambda_val)> lambda_vec(lambda_val);
  const size_t N = max_size(y, mu, sigma, lambda);
  T_partials_return log_cdf = 0.0;

  for (size_t n = 0; n < N; ++n) {
    if (value_of_rec(y_vec[n]) == INFTY) {
      continue;
    }
    const T_partials_return sigma_dbl = sigma_vec[n];
    const T_partials_return inv_sigma = 1.0 / sigma_dbl;
    const T_partials_return z = (y_vec[n] - mu_vec[n]) * inv_sigma;
    const T_partials_return a = lambda_vec[n] * sigma_dbl;
    const auto terms = internal::exp_mod_normal_cdf_terms(z, a);
    log_cdf += terms.log_cdf;

    if constexpr (is_any_autodiff_v<T_y, T_loc, T_scale, T_inv_scale>) {
      if constexpr (is_autodiff_v<T_y>) {
        partials<0>(ops_partials)[n] += terms.dz_log_cdf * inv_sigma;
      }
      if constexpr (is_autodiff_v<T_loc>) {
        partials<1>(ops_partials)[n] -= terms.dz_log_cdf * inv_sigma;
      }
      if constexpr (is_autodiff_v<T_scale>) {
        partials<2>(ops_partials)[n]
            += (-z * terms.dz_log_cdf + terms.a_da_log_cdf) * inv_sigma;
      }
    }
  }

  using std::exp;
  const T_partials_return cdf = exp(log_cdf);
  if constexpr (is_autodiff_v<T_y>) {
    for (size_t n = 0; n < stan::math::size(y); ++n) {
      partials<0>(ops_partials)[n] *= cdf;
    }
  }
  if constexpr (is_autodiff_v<T_loc>) {
    for (size_t n = 0; n < stan::math::size(mu); ++n) {
      partials<1>(ops_partials)[n] *= cdf;
    }
  }
  if constexpr (is_autodiff_v<T_scale>) {
    for (size_t n = 0; n < stan::math::size(sigma); ++n) {
      partials<2>(ops_partials)[n] *= cdf;
    }
  }
  if constexpr (is_autodiff_v<T_inv_scale>) {
    for (size_t n = 0; n < N; ++n) {
      if (value_of_rec(y_vec[n]) == INFTY) {
        continue;
      }
      const T_partials_return sigma_dbl = sigma_vec[n];
      const T_partials_return z = (y_vec[n] - mu_vec[n]) / sigma_dbl;
      const T_partials_return a = lambda_vec[n] * sigma_dbl;
      const auto terms = internal::exp_mod_normal_cdf_terms(z, a);
      partials<3>(ops_partials)[n]
          += exp(log_cdf - log(a)) * sigma_dbl * terms.a_da_log_cdf;
    }
  }
  return ops_partials.build(cdf);
}

}  // namespace math
}  // namespace stan
#endif
