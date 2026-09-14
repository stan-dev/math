#ifndef STAN_MATH_PRIM_PROB_EXP_MOD_NORMAL_LCCDF_HPP
#define STAN_MATH_PRIM_PROB_EXP_MOD_NORMAL_LCCDF_HPP

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
inline return_type_t<T_y, T_loc, T_scale, T_inv_scale> exp_mod_normal_lccdf(
    const T_y& y, const T_loc& mu, const T_scale& sigma,
    const T_inv_scale& lambda) {
  using T_partials_return = partials_return_t<T_y, T_loc, T_scale, T_inv_scale>;
  using T_y_ref = ref_type_if_not_constant_t<T_y>;
  using T_mu_ref = ref_type_if_not_constant_t<T_loc>;
  using T_sigma_ref = ref_type_if_not_constant_t<T_scale>;
  using T_lambda_ref = ref_type_if_not_constant_t<T_inv_scale>;
  static constexpr const char* function = "exp_mod_normal_lccdf";
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
    return 0;
  }

  auto ops_partials
      = make_partials_propagator(y_ref, mu_ref, sigma_ref, lambda_ref);

  scalar_seq_view<decltype(y_val)> y_vec(y_val);
  scalar_seq_view<decltype(mu_val)> mu_vec(mu_val);
  scalar_seq_view<decltype(sigma_val)> sigma_vec(sigma_val);
  scalar_seq_view<decltype(lambda_val)> lambda_vec(lambda_val);
  const size_t N = max_size(y, mu, sigma, lambda);

  for (size_t n = 0, size_y = stan::math::size(y); n < size_y; ++n) {
    if (value_of_rec(y_vec[n]) == INFTY) {
      return ops_partials.build(negative_infinity());
    }
  }

  T_partials_return ccdf_log(0.0);
  for (size_t n = 0; n < N; ++n) {
    if (value_of_rec(y_vec[n]) == NEGATIVE_INFTY) {
      continue;
    }
    const T_partials_return sigma_dbl = sigma_vec[n];
    const T_partials_return inv_sigma = 1.0 / sigma_dbl;
    const T_partials_return z = (y_vec[n] - mu_vec[n]) * inv_sigma;
    const T_partials_return a = lambda_vec[n] * sigma_dbl;
    const auto terms = internal::exp_mod_normal_cdf_terms(z, a);
    ccdf_log += terms.log_ccdf;

    if constexpr (is_any_autodiff_v<T_y, T_loc, T_scale, T_inv_scale>) {
      if constexpr (is_autodiff_v<T_y>) {
        partials<0>(ops_partials)[n] += terms.dz_log_ccdf * inv_sigma;
      }
      if constexpr (is_autodiff_v<T_loc>) {
        partials<1>(ops_partials)[n] -= terms.dz_log_ccdf * inv_sigma;
      }
      if constexpr (is_autodiff_v<T_scale>) {
        partials<2>(ops_partials)[n]
            += (-z * terms.dz_log_ccdf + a * terms.da_log_ccdf) * inv_sigma;
      }
      if constexpr (is_autodiff_v<T_inv_scale>) {
        partials<3>(ops_partials)[n] += sigma_dbl * terms.da_log_ccdf;
      }
    }
  }

  return ops_partials.build(ccdf_log);
}

}  // namespace math
}  // namespace stan
#endif
