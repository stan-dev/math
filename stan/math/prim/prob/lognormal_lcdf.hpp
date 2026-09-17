#ifndef STAN_MATH_PRIM_PROB_LOGNORMAL_LCDF_HPP
#define STAN_MATH_PRIM_PROB_LOGNORMAL_LCDF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/as_column_vector_or_scalar.hpp>
#include <stan/math/prim/fun/as_array_or_scalar.hpp>
#include <stan/math/prim/fun/as_value_column_array_or_scalar.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/log.hpp>
#include <stan/math/prim/fun/max_size.hpp>
#include <stan/math/prim/fun/promote_scalar.hpp>
#include <stan/math/prim/fun/size.hpp>
#include <stan/math/prim/fun/size_zero.hpp>
#include <stan/math/prim/fun/to_ref.hpp>
#include <stan/math/prim/fun/value_of.hpp>
#include <stan/math/prim/functor/apply_scalar_binary.hpp>
#include <stan/math/prim/functor/partials_propagator.hpp>
#include <stan/math/prim/prob/std_normal_lcdf_impl.hpp>
#include <cmath>

namespace stan {
namespace math {

template <typename T_y, typename T_loc, typename T_scale,
          require_all_not_nonscalar_prim_or_rev_kernel_expression_t<
              T_y, T_loc, T_scale>* = nullptr>
inline return_type_t<T_y, T_loc, T_scale> lognormal_lcdf(const T_y& y,
                                                         const T_loc& mu,
                                                         const T_scale& sigma) {
  using T_partials_return = partials_return_t<T_y, T_loc, T_scale>;
  using T_y_ref = ref_type_if_not_constant_t<T_y>;
  using T_mu_ref = ref_type_if_not_constant_t<T_loc>;
  using T_sigma_ref = ref_type_if_not_constant_t<T_scale>;
  static constexpr const char* function = "lognormal_lcdf";

  check_consistent_sizes(function, "Random variable", y, "Location parameter",
                         mu, "Scale parameter", sigma);
  T_y_ref y_ref = y;
  T_mu_ref mu_ref = mu;
  T_sigma_ref sigma_ref = sigma;

  decltype(auto) y_val = to_ref(as_value_column_array_or_scalar(y_ref));
  decltype(auto) mu_val = to_ref(as_value_column_array_or_scalar(mu_ref));
  decltype(auto) sigma_val = to_ref(as_value_column_array_or_scalar(sigma_ref));

  check_nonnegative(function, "Random variable", y_val);
  check_finite(function, "Location parameter", mu_val);
  check_positive_finite(function, "Scale parameter", sigma_val);

  if (size_zero(y, mu, sigma)) {
    return 0;
  }

  auto ops_partials = make_partials_propagator(y_ref, mu_ref, sigma_ref);

  if (sum(promote_scalar<int>(y_val == 0))) {
    return ops_partials.build(NEGATIVE_INFTY);
  }

  const auto& z = to_ref((log(y_val) - mu_val) / sigma_val);
  const auto [values, slopes] = internal::std_normal_lcdf_value_grad<
      is_any_autodiff_v<T_y, T_loc, T_scale>>(z);
  T_partials_return cdf_log = sum(values);

  if constexpr (is_any_autodiff_v<T_y, T_loc, T_scale>) {
    const auto& scaled_slope = to_ref(slopes / sigma_val);
    if constexpr (is_autodiff_v<T_y>) {
      partials<0>(ops_partials) = scaled_slope / y_val;
    }
    if constexpr (is_autodiff_v<T_loc>) {
      partials<1>(ops_partials) = -scaled_slope;
    }
    if constexpr (is_autodiff_v<T_scale>) {
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
