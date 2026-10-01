#ifndef STAN_MATH_PRIM_PROB_SKEW_NORMAL_LPDF_HPP
#define STAN_MATH_PRIM_PROB_SKEW_NORMAL_LPDF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/as_column_vector_or_scalar.hpp>
#include <stan/math/prim/fun/as_array_or_scalar.hpp>
#include <stan/math/prim/fun/as_value_column_array_or_scalar.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/inv.hpp>
#include <stan/math/prim/fun/log.hpp>
#include <stan/math/prim/fun/max_size.hpp>
#include <stan/math/prim/fun/size.hpp>
#include <stan/math/prim/fun/size_zero.hpp>
#include <stan/math/prim/fun/square.hpp>
#include <stan/math/prim/fun/to_ref.hpp>
#include <stan/math/prim/fun/value_of.hpp>
#include <stan/math/prim/functor/partials_propagator.hpp>
#include <stan/math/prim/fun/std_normal_lcdf_impl.hpp>
#include <cmath>

namespace stan {
namespace math {

template <bool propto, typename T_y, typename T_loc, typename T_scale,
          typename T_shape,
          require_all_not_nonscalar_prim_or_rev_kernel_expression_t<
              T_y, T_loc, T_scale, T_shape>* = nullptr>
inline return_type_t<T_y, T_loc, T_scale, T_shape> skew_normal_lpdf(
    T_y&& y, T_loc&& mu, T_scale&& sigma, T_shape&& alpha) {
  using T_partials_return = partials_return_t<T_y, T_loc, T_scale, T_shape>;
  using T_y_ref = ref_type_if_not_constant_t<T_y>;
  using T_mu_ref = ref_type_if_not_constant_t<T_loc>;
  using T_sigma_ref = ref_type_if_not_constant_t<T_scale>;
  using T_alpha_ref = ref_type_if_not_constant_t<T_shape>;
  static constexpr const char* function = "skew_normal_lpdf";
  check_consistent_sizes(function, "Random variable", y, "Location parameter",
                         mu, "Scale parameter", sigma, "Shape parameter",
                         alpha);
  T_y_ref y_ref = std::forward<T_y>(y);
  T_mu_ref mu_ref = std::forward<T_loc>(mu);
  T_sigma_ref sigma_ref = std::forward<T_scale>(sigma);
  T_alpha_ref alpha_ref = std::forward<T_shape>(alpha);

  decltype(auto) y_val = to_ref(as_value_column_array_or_scalar(y_ref));
  decltype(auto) mu_val = to_ref(as_value_column_array_or_scalar(mu_ref));
  decltype(auto) sigma_val = to_ref(as_value_column_array_or_scalar(sigma_ref));
  decltype(auto) alpha_val = to_ref(as_value_column_array_or_scalar(alpha_ref));

  check_not_nan(function, "Random variable", y_val);
  check_finite(function, "Location parameter", mu_val);
  check_finite(function, "Shape parameter", alpha_val);
  check_positive(function, "Scale parameter", sigma_val);

  if (size_zero(y_ref, mu_ref, sigma_ref, alpha_ref)) {
    return 0.0;
  }
  if constexpr (!include_summand<propto, T_y, T_loc, T_scale, T_shape>::value) {
    return 0.0;
  }

  auto ops_partials
      = make_partials_propagator(y_ref, mu_ref, sigma_ref, alpha_ref);

  const auto& inv_sigma
      = to_ref_if<is_any_autodiff_v<T_y, T_loc, T_scale>>(inv(sigma_val));
  const auto& z = to_ref((y_val - mu_val) * inv_sigma);
  const auto& az = to_ref(alpha_val * z);
  const auto [values, slopes] = internal::std_normal_lcdf_value_grad<
      is_any_autodiff_v<T_y, T_loc, T_scale, T_shape>>(az);

  size_t N = max_size(y_ref, mu_ref, sigma_ref, alpha_ref);
  T_partials_return logp = N * LOG_TWO + sum(values);
  if constexpr (include_summand<propto>::value) {
    logp -= HALF_LOG_TWO_PI * N;
  }
  if constexpr (include_summand<propto, T_scale>::value) {
    logp -= sum(log(sigma_val)) * N / math::size(sigma_ref);
  }
  if constexpr (include_summand<propto, T_y, T_loc, T_scale>::value) {
    logp -= sum(square(z)) * 0.5 * N / max_size(y_ref, mu_ref, sigma_ref);
  }
  if constexpr (is_any_autodiff_v<T_y, T_loc, T_scale>) {
    const auto& score = to_ref_if<
        (is_autodiff_v<T_y> + is_autodiff_v<T_loc> + is_autodiff_v<T_scale>)
        >= 2>((slopes * alpha_val - z) * inv_sigma);
    if constexpr (is_autodiff_v<T_y>) {
      partials<0>(ops_partials) = score;
    }
    if constexpr (is_autodiff_v<T_loc>) {
      partials<1>(ops_partials) = -score;
    }
    if constexpr (is_autodiff_v<T_scale>) {
      partials<2>(ops_partials) = -score * z - inv_sigma;
    }
  }
  if constexpr (is_autodiff_v<T_shape>) {
    partials<3>(ops_partials) = slopes * z;
  }
  return ops_partials.build(logp);
}

template <typename T_y, typename T_loc, typename T_scale, typename T_shape>
inline return_type_t<T_y, T_loc, T_scale, T_shape> skew_normal_lpdf(
    T_y&& y, T_loc&& mu, T_scale&& sigma, T_shape&& alpha) {
  return skew_normal_lpdf<false>(std::forward<T_y>(y), std::forward<T_loc>(mu),
                                 std::forward<T_scale>(sigma),
                                 std::forward<T_shape>(alpha));
}

}  // namespace math
}  // namespace stan
#endif
