#ifndef STAN_MATH_PRIM_PROB_VON_MISES_LPDF_HPP
#define STAN_MATH_PRIM_PROB_VON_MISES_LPDF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/as_value_column_array_or_scalar.hpp>
#include <stan/math/prim/fun/bessel_i0e_derivative_ratio.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/log_bessel_i0e.hpp>
#include <stan/math/prim/fun/max_size.hpp>
#include <stan/math/prim/fun/sin.hpp>
#include <stan/math/prim/fun/size.hpp>
#include <stan/math/prim/fun/size_zero.hpp>
#include <stan/math/prim/fun/square.hpp>
#include <stan/math/prim/fun/to_ref.hpp>
#include <stan/math/prim/functor/partials_propagator.hpp>
#include <cmath>
#include <utility>

namespace stan {
namespace math {

template <bool propto, typename T_y, typename T_loc, typename T_scale>
inline return_type_t<T_y, T_loc, T_scale> von_mises_lpdf(T_y const& y,
                                                         T_loc const& mu,
                                                         T_scale const& kappa) {
  using T_partials_return = partials_return_t<T_y, T_loc, T_scale>;
  using T_y_ref = ref_type_if_not_constant_t<T_y>;
  using T_mu_ref = ref_type_if_not_constant_t<T_loc>;
  using T_kappa_ref = ref_type_if_not_constant_t<T_scale>;
  static char const* const function = "von_mises_lpdf";
  check_consistent_sizes(function, "Random variable", y, "Location parameter",
                         mu, "Scale parameter", kappa);
  T_y_ref y_ref = y;
  T_mu_ref mu_ref = mu;
  T_kappa_ref kappa_ref = kappa;

  decltype(auto) y_val = to_ref(as_value_column_array_or_scalar(y_ref));
  decltype(auto) mu_val = to_ref(as_value_column_array_or_scalar(mu_ref));
  decltype(auto) kappa_val = to_ref(as_value_column_array_or_scalar(kappa_ref));
  check_finite(function, "Random variable", y_val);
  check_finite(function, "Location parameter", mu_val);
  check_nonnegative(function, "Scale parameter", kappa_val);
  check_finite(function, "Scale parameter", kappa_val);

  if (size_zero(y, mu, kappa)) {
    return 0;
  }
  if constexpr (!include_summand<propto, T_y, T_loc, T_scale>::value) {
    return 0;
  }

  auto ops_partials = make_partials_propagator(y_ref, mu_ref, kappa_ref);

  // 1 - cos(y - mu), without cancellation for y close to mu
  const auto& one_minus_cos = to_ref(2.0 * square(sin(0.5 * (y_val - mu_val))));

  size_t N = max_size(y, mu, kappa);
  T_partials_return logp = 0.0 - sum(kappa_val * one_minus_cos);
  if constexpr (include_summand<propto>::value) {
    logp -= LOG_TWO_PI * N;
  }
  if constexpr (include_summand<propto, T_scale>::value) {
    logp -= sum(log_bessel_i0e(kappa_val)) * N / math::size(kappa);
  }
  if constexpr (is_autodiff_v<T_scale>) {
    partials<2>(ops_partials)
        = bessel_i0e_derivative_ratio(0, kappa_val) - one_minus_cos;
  }

  if constexpr (is_any_autodiff_v<T_y, T_loc>) {
    const auto& sin_diff = sin(y_val - mu_val);
    auto kappa_sin = to_ref_if<(is_autodiff_v<T_y> && is_autodiff_v<T_loc>)>(
        kappa_val * sin_diff);
    if constexpr (is_autodiff_v<T_y>) {
      partials<0>(ops_partials) = -kappa_sin;
    }
    if constexpr (is_autodiff_v<T_loc>) {
      partials<1>(ops_partials) = std::move(kappa_sin);
    }
  }

  return ops_partials.build(logp);
}

template <typename T_y, typename T_loc, typename T_scale>
inline return_type_t<T_y, T_loc, T_scale> von_mises_lpdf(T_y const& y,
                                                         T_loc const& mu,
                                                         T_scale const& kappa) {
  return von_mises_lpdf<false>(y, mu, kappa);
}

}  // namespace math
}  // namespace stan
#endif
