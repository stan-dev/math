#ifndef STAN_MATH_PRIM_PROB_NORMAL_LCDF_HPP
#define STAN_MATH_PRIM_PROB_NORMAL_LCDF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/as_value_column_array_or_scalar.hpp>
#include <stan/math/prim/fun/sum.hpp>
#include <stan/math/prim/fun/to_ref.hpp>
#include <stan/math/prim/fun/select.hpp>
#include <stan/math/prim/fun/size_zero.hpp>
#include <stan/math/prim/functor/partials_propagator.hpp>
#include <stan/math/prim/prob/std_normal_lcdf_impl.hpp>

namespace stan {
namespace math {
namespace internal {

/** Log of the normal cdf, or of its complement when `reflect` is set: the
 * standardized value is negated so no negated autodiff operands are built.
 */
template <bool reflect, typename T_y, typename T_loc, typename T_scale>
inline return_type_t<T_y, T_loc, T_scale> normal_lcdf_impl(const char* function,
                                                           T_y&& y, T_loc&& mu,
                                                           T_scale&& sigma) {
  using T_partials_return = partials_return_t<T_y, T_loc, T_scale>;
  using T_y_ref = ref_type_if_not_constant_t<T_y>;
  using T_mu_ref = ref_type_if_not_constant_t<T_loc>;
  using T_sigma_ref = ref_type_if_not_constant_t<T_scale>;
  constexpr double sign = reflect ? -1.0 : 1.0;
  check_consistent_sizes(function, "Random variable", y, "Location parameter",
                         mu, "Scale parameter", sigma);
  T_y_ref y_ref = std::forward<T_y>(y);
  T_mu_ref mu_ref = std::forward<T_loc>(mu);
  T_sigma_ref sigma_ref = std::forward<T_scale>(sigma);
  decltype(auto) y_val = to_ref(as_value_column_array_or_scalar(y_ref));
  decltype(auto) mu_val = to_ref(as_value_column_array_or_scalar(mu_ref));
  decltype(auto) sigma_val = to_ref(as_value_column_array_or_scalar(sigma_ref));
  check_not_nan(function, "Random variable", y_val);
  check_finite(function, "Location parameter", mu_val);
  check_positive(function, "Scale parameter", sigma_val);

  if (size_zero(y_ref, mu_ref, sigma_ref)) {
    return 0;
  }

  const auto& z = to_ref(sign * (y_val - mu_val) / sigma_val);
  const auto [values, slopes] = internal::std_normal_lcdf_value_grad<
      is_any_autodiff_v<T_y, T_loc, T_scale>>(z);
  const T_partials_return cdf_log = sum(values);
  auto ops_partials = make_partials_propagator(y_ref, mu_ref, sigma_ref);
  if constexpr (is_any_autodiff_v<T_y, T_loc, T_scale>) {
    const auto& scaled_slope
        = to_ref_if<(is_autodiff_v<T_y> + is_autodiff_v<T_loc>
                     + is_autodiff_v<T_scale>)
                    >= 2>(slopes / sigma_val);
    if constexpr (is_autodiff_v<T_y>) {
      partials<0>(ops_partials) = sign * scaled_slope;
    }
    if constexpr (is_autodiff_v<T_loc>) {
      partials<1>(ops_partials) = -sign * scaled_slope;
    }
    if constexpr (is_autodiff_v<T_scale>) {
      // The positive infinite endpoint has a zero slope, not 0 * infinity.
      partials<2>(ops_partials)
          = select(scaled_slope == 0, 0.0, -scaled_slope * z);
    }
  }
  return ops_partials.build(cdf_log);
}

}  // namespace internal

/** \ingroup prob_dists
 * @brief Calculates the log of the cdf of the normal distribution
 *
 * Uses the shared standard-normal kernel after standardization. See
 * std_normal_lcdf_impl.hpp for the Cody approximations and their crossovers.
 *
 * @tparam T_y A vector or scalar type for the random variable.
 * @tparam T_loc A vector or scalar type for the location parameter.
 * @tparam T_scale A vector or scalar type for the scale parameter.
 * @param y (Sequence of) scalar(s).
 * @param mu (Sequence of) scalar(s).
 * @param sigma (Sequence of) scalar(s).
 * @return The log of the normal cdf evaluated at the specified arguments. If
 *   given containers, the log of the product of the cdfs.
 */
template <typename T_y, typename T_loc, typename T_scale,
          require_all_not_nonscalar_prim_or_rev_kernel_expression_t<
              T_y, T_loc, T_scale>* = nullptr>
inline return_type_t<T_y, T_loc, T_scale> normal_lcdf(T_y&& y, T_loc&& mu,
                                                      T_scale&& sigma) {
  return internal::normal_lcdf_impl<false>("normal_lcdf", std::forward<T_y>(y),
                                           std::forward<T_loc>(mu),
                                           std::forward<T_scale>(sigma));
}

}  // namespace math
}  // namespace stan
#endif
