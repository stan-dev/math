#ifndef STAN_MATH_PRIM_PROB_NORMAL_CDF_HPP
#define STAN_MATH_PRIM_PROB_NORMAL_CDF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/as_value_column_array_or_scalar.hpp>
#include <stan/math/prim/fun/eval.hpp>
#include <stan/math/prim/fun/prod.hpp>
#include <stan/math/prim/fun/select.hpp>
#include <stan/math/prim/fun/size_zero.hpp>
#include <stan/math/prim/fun/std_normal_lcdf_impl.hpp>
#include <stan/math/prim/fun/to_ref.hpp>
#include <stan/math/prim/functor/partials_propagator.hpp>

namespace stan {
namespace math {

/** \ingroup prob_dists
 * Calculates the normal cumulative distribution function for the given
 * variate, location, and scale.
 *
 * \f$\Phi(x) = \frac{1}{\sqrt{2 \pi}} \int_{-\inf}^x e^{-t^2/2} dt\f$.
 *
 * Each factor is `0.5 * erfc(-z / sqrt(2))` with `z = (y - mu) / sigma`,
 * which has no cancellation for either sign; the partials are the cdf times
 * `phi(z) / Phi(z)`.
 *
 * @tparam T_y type of y
 * @tparam T_loc type of mean parameter
 * @tparam T_scale type of standard deviation parameter
 * @param y A scalar variate.
 * @param mu The location of the normal distribution.
 * @param sigma The scale of the normal distribution
 * @return The unit normal cdf evaluated at the specified arguments.
 */
template <typename T_y, typename T_loc, typename T_scale,
          require_all_not_nonscalar_prim_or_rev_kernel_expression_t<
              T_y, T_loc, T_scale>* = nullptr>
inline return_type_t<T_y, T_loc, T_scale> normal_cdf(T_y&& y, T_loc&& mu,
                                                     T_scale&& sigma) {
  using T_y_ref = ref_type_if_not_constant_t<T_y>;
  using T_mu_ref = ref_type_if_not_constant_t<T_loc>;
  using T_sigma_ref = ref_type_if_not_constant_t<T_scale>;
  static constexpr const char* function = "normal_cdf";
  check_consistent_sizes(function, "Random variable", y, "Location parameter",
                         mu, "Scale parameter", sigma);
  T_y_ref y_ref = std::forward<T_y>(y);
  T_mu_ref mu_ref = std::forward<T_loc>(mu);
  T_sigma_ref sigma_ref = std::forward<T_scale>(sigma);
  // to_ref only reduces compile memory by evaluating expression arguments.
  decltype(auto) y_val = to_ref(as_value_column_array_or_scalar(y_ref));
  decltype(auto) mu_val = to_ref(as_value_column_array_or_scalar(mu_ref));
  decltype(auto) sigma_val = to_ref(as_value_column_array_or_scalar(sigma_ref));
  check_not_nan(function, "Random variable", y_val);
  check_finite(function, "Location parameter", mu_val);
  check_positive(function, "Scale parameter", sigma_val);

  if (size_zero(y_ref, mu_ref, sigma_ref)) {
    return 1.0;
  }

  auto z = (y_val - mu_val) / sigma_val;
  if constexpr (is_any_autodiff_v<T_y, T_loc, T_scale>) {
    const auto [values, ratios] = internal::std_normal_cdf_value_ratio<true>(z);
    const auto cdf = prod(values);
    auto ops_partials = make_partials_propagator(y_ref, mu_ref, sigma_ref);
    auto scaled = eval(cdf * ratios / sigma_val);
    if constexpr (is_autodiff_v<T_y>) {
      partials<0>(ops_partials) = scaled;
    }
    if constexpr (is_autodiff_v<T_loc>) {
      partials<1>(ops_partials) = -scaled;
    }
    if constexpr (is_autodiff_v<T_scale>) {
      // The infinite endpoints have a zero partial, not 0 * infinity.
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
