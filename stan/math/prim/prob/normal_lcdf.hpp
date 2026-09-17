#ifndef STAN_MATH_PRIM_PROB_NORMAL_LCDF_HPP
#define STAN_MATH_PRIM_PROB_NORMAL_LCDF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/max_size.hpp>
#include <stan/math/prim/fun/as_value_column_array_or_scalar.hpp>
#include <stan/math/prim/fun/sum.hpp>
#include <stan/math/prim/fun/to_ref.hpp>
#include <stan/math/prim/functor/apply_scalar_binary.hpp>
#include <stan/math/prim/fun/size_zero.hpp>
#include <stan/math/prim/functor/partials_propagator.hpp>
#include <stan/math/prim/prob/std_normal_lcdf_impl.hpp>
#include <stan/math/prim/prob/normal_standardize.hpp>

namespace stan {
namespace math {
namespace internal {
constexpr char normal_lcdf_func[] = "normal_lcdf";
}  // namespace internal

/** \ingroup prob_dists
 * @brief Calculates the log of the cdf of the normal distribution
 *
 * Uses the shared standard-normal tail kernel after standardization.
 * See std_normal_lcdf_impl.hpp for the Cody approximation and its crossover.
 * Values and slopes share their erfc or tail evaluation.
 *
 * @tparam func name reported by the error checks. Reflected distributions
 *   such as `normal_lccdf` delegate here and pass their own name so that
 *   exceptions name the function the user actually called.
 * @tparam T_y A vector or scalar type for the random variable.
 * @tparam T_loc A vector or scalar type for the location parameter.
 * @tparam T_scale A vector or scalar type for the scale parameter.
 * @param y (Sequence of) scalar(s).
 * @param mu (Sequence of) scalar(s).
 * @param sigma (Sequence of) scalar(s).
 * @return The log of the normal cdf evaluated at the specified arguments. If
 *   given containers, the log of the product of the cdfs.
 */
template <const char* func = internal::normal_lcdf_func, typename T_y,
          typename T_loc, typename T_scale,
          require_all_not_nonscalar_prim_or_rev_kernel_expression_t<
              T_y, T_loc, T_scale>* = nullptr>
inline return_type_t<T_y, T_loc, T_scale> normal_lcdf(const T_y& y,
                                                      const T_loc& mu,
                                                      const T_scale& sigma) {
  using T_partials_return = partials_return_t<T_y, T_loc, T_scale>;
  using T_y_ref = ref_type_t<T_y>;
  using T_mu_ref = ref_type_t<T_loc>;
  using T_sigma_ref = ref_type_t<T_scale>;
  static constexpr const char* function = func;
  check_consistent_sizes(function, "Random variable", y, "Location parameter",
                         mu, "Scale parameter", sigma);
  T_y_ref y_ref = y;
  T_mu_ref mu_ref = mu;
  T_sigma_ref sigma_ref = sigma;
  check_not_nan(function, "Random variable", y_ref);
  check_finite(function, "Location parameter", mu_ref);
  check_positive(function, "Scale parameter", sigma_ref);

  if (size_zero(y, mu, sigma)) {
    return 0;
  }

  const auto& y_val = to_ref(as_value_column_array_or_scalar(y_ref));
  const auto& mu_val = to_ref(as_value_column_array_or_scalar(mu_ref));
  const auto& sigma_val = to_ref(as_value_column_array_or_scalar(sigma_ref));
  const auto& z
      = to_ref(internal::normal_standardize(y_val, mu_val, sigma_val));
  const auto [values, slopes] = internal::std_normal_lcdf_value_grad<
      is_any_autodiff_v<T_y, T_loc, T_scale>>(z);
  const T_partials_return cdf_log = sum(values);
  auto ops_partials = make_partials_propagator(y_ref, mu_ref, sigma_ref);
  if constexpr (is_any_autodiff_v<T_y, T_loc, T_scale>) {
    const auto& scaled_slope = to_ref(slopes / sigma_val);
    if constexpr (is_autodiff_v<T_y>) {
      partials<0>(ops_partials) = scaled_slope;
    }
    if constexpr (is_autodiff_v<T_loc>) {
      partials<1>(ops_partials) = -scaled_slope;
    }
    if constexpr (is_autodiff_v<T_scale>) {
      // The positive infinite endpoint has a zero slope, not 0 * infinity.
      partials<2>(ops_partials) = apply_scalar_binary(
          [](const auto& slope, const auto& zi) {
            using R = return_type_t<decltype(slope), decltype(zi)>;
            return slope == 0 ? R(0) : R(-slope * zi);
          },
          scaled_slope, z);
    }
  }
  return ops_partials.build(cdf_log);
}

}  // namespace math
}  // namespace stan
#endif
