#ifndef STAN_MATH_PRIM_PROB_LOGNORMAL_CDF_HPP
#define STAN_MATH_PRIM_PROB_LOGNORMAL_CDF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/any.hpp>
#include <stan/math/prim/fun/as_value_column_array_or_scalar.hpp>
#include <stan/math/prim/fun/eval.hpp>
#include <stan/math/prim/fun/log.hpp>
#include <stan/math/prim/fun/prod.hpp>
#include <stan/math/prim/fun/select.hpp>
#include <stan/math/prim/fun/size_zero.hpp>
#include <stan/math/prim/fun/std_normal_lcdf_impl.hpp>
#include <stan/math/prim/fun/to_ref.hpp>
#include <stan/math/prim/functor/partials_propagator.hpp>

namespace stan {
namespace math {

/** \ingroup prob_dists
 * The lognormal cumulative distribution function. Each factor is
 * `0.5 * erfc(-z / sqrt(2))` with `z = (log(y) - mu) / sigma`, which has no
 * cancellation for either sign; the partials are the cdf times
 * `phi(z) / Phi(z)`.
 *
 * @tparam T_y type of y
 * @tparam T_loc type of location parameter
 * @tparam T_scale type of scale parameter
 * @param y A scalar variate.
 * @param mu The location of the lognormal distribution.
 * @param sigma The scale of the lognormal distribution.
 * @return The lognormal cdf evaluated at the specified arguments.
 */
template <typename T_y, typename T_loc, typename T_scale,
          require_all_not_nonscalar_prim_or_rev_kernel_expression_t<
              T_y, T_loc, T_scale>* = nullptr>
inline return_type_t<T_y, T_loc, T_scale> lognormal_cdf(T_y&& y, T_loc&& mu,
                                                        T_scale&& sigma) {
  using T_y_ref = ref_type_if_not_constant_t<T_y>;
  using T_mu_ref = ref_type_if_not_constant_t<T_loc>;
  using T_sigma_ref = ref_type_if_not_constant_t<T_scale>;
  static constexpr const char* function = "lognormal_cdf";
  check_consistent_sizes(function, "Random variable", y, "Location parameter",
                         mu, "Scale parameter", sigma);
  T_y_ref y_ref = std::forward<T_y>(y);
  T_mu_ref mu_ref = std::forward<T_loc>(mu);
  T_sigma_ref sigma_ref = std::forward<T_scale>(sigma);
  // to_ref only reduces compile memory by evaluating expression arguments.
  decltype(auto) y_val = to_ref(as_value_column_array_or_scalar(y_ref));
  decltype(auto) mu_val = to_ref(as_value_column_array_or_scalar(mu_ref));
  decltype(auto) sigma_val = to_ref(as_value_column_array_or_scalar(sigma_ref));
  check_nonnegative(function, "Random variable", y_val);
  check_finite(function, "Location parameter", mu_val);
  check_positive_finite(function, "Scale parameter", sigma_val);

  if (size_zero(y_ref, mu_ref, sigma_ref)) {
    return 1.0;
  }

  auto ops_partials = make_partials_propagator(y_ref, mu_ref, sigma_ref);
  if (any(y_val == 0)) {
    return ops_partials.build(0.0);
  }

  auto z = (log(y_val) - mu_val) / sigma_val;
  if constexpr (is_any_autodiff_v<T_y, T_loc, T_scale>) {
    const auto [values, ratios] = internal::std_normal_cdf_value_ratio<true>(z);
    const auto cdf = prod(values);
    auto scaled = eval(cdf * ratios / sigma_val);
    if constexpr (is_autodiff_v<T_y>) {
      partials<0>(ops_partials) = scaled / y_val;
    }
    if constexpr (is_autodiff_v<T_loc>) {
      partials<1>(ops_partials) = -scaled;
    }
    if constexpr (is_autodiff_v<T_scale>) {
      // The infinite endpoint has a zero partial, not 0 * infinity.
      partials<2>(ops_partials) = select(scaled == 0, 0.0, -scaled * z);
    }
    return ops_partials.build(cdf);
  } else {
    return prod(internal::std_normal_cdf_value_ratio<false>(std::move(z)));
  }
}

}  // namespace math
}  // namespace stan
#endif
