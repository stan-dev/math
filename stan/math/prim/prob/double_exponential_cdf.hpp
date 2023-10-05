#ifndef STAN_MATH_PRIM_PROB_DOUBLE_EXPONENTIAL_CDF_HPP
#define STAN_MATH_PRIM_PROB_DOUBLE_EXPONENTIAL_CDF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/as_column_vector_or_scalar.hpp>
#include <stan/math/prim/fun/as_array_or_scalar.hpp>
#include <stan/math/prim/fun/as_value_column_array_or_scalar.hpp>
#include <stan/math/prim/fun/eigen_comparisons.hpp>
#include <stan/math/prim/fun/exp.hpp>
#include <stan/math/prim/fun/inv.hpp>
#include <stan/math/prim/fun/max_size.hpp>
#include <stan/math/prim/fun/size_zero.hpp>
#include <stan/math/prim/fun/to_ref.hpp>
#include <stan/math/prim/fun/value_of.hpp>
#include <stan/math/prim/functor/partials_propagator.hpp>
#include <cmath>

namespace stan {
namespace math {

/** \ingroup prob_dists
 * Returns the double exponential cumulative density function. Given
 * containers of matching sizes, returns the product of probabilities.
 *
 * @tparam T_y type of real parameter.
 * @tparam T_loc type of location parameter.
 * @tparam T_scale type of scale parameter.
 * @param y real parameter
 * @param mu location parameter
 * @param sigma scale parameter
 * @return probability or product of probabilities
 * @throw std::domain_error if y is nan, mu is infinite,
 *  or sigma is nonpositive
 */
template <typename T_y, typename T_loc, typename T_scale,
          require_all_not_nonscalar_prim_or_rev_kernel_expression_t<
              T_y, T_loc, T_scale>* = nullptr>
return_type_t<T_y, T_loc, T_scale> double_exponential_cdf(
    const T_y& y, const T_loc& mu, const T_scale& sigma) {
  using T_partials_return = partials_return_t<T_y, T_loc, T_scale>;
  using T_partials_array = Eigen::Array<T_partials_return, Eigen::Dynamic, 1>;
  using T_rep_deriv
      = std::conditional_t<is_vector<T_y>::value || is_vector<T_loc>::value
                               || is_vector<T_scale>::value,
                           T_partials_array, T_partials_return>;
  using std::exp;
  using T_y_ref = ref_type_if_not_constant_t<T_y>;
  using T_mu_ref = ref_type_if_not_constant_t<T_loc>;
  using T_sigma_ref = ref_type_if_not_constant_t<T_scale>;
  static const char* function = "double_exponential_cdf";
  T_y_ref y_ref = y;
  T_mu_ref mu_ref = mu;
  T_sigma_ref sigma_ref = sigma;

  T_partials_return cdf(1.0);
  auto ops_partials = make_partials_propagator(y_ref, mu_ref, sigma_ref);

  decltype(auto) y_val = to_ref(as_value_column_array_or_scalar(y_ref));
  decltype(auto) mu_val = to_ref(as_value_column_array_or_scalar(mu_ref));
  decltype(auto) sigma_val = to_ref(as_value_column_array_or_scalar(sigma_ref));

  check_not_nan(function, "Random variable", y_val);
  check_finite(function, "Location parameter", mu_val);
  check_positive_finite(function, "Scale parameter", sigma_val);

  if (size_zero(y, mu, sigma)) {
    return 1.0;
  }

  const auto& inv_sigma = to_ref(inv(sigma_val));
  const auto& scaled_diff = to_ref_if<!is_constant_all<T_scale>::value>(
      (y_val - mu_val) * inv_sigma);
  const auto& exp_scaled_diff = to_ref(exp(scaled_diff));

  T_rep_deriv rep_deriv;
  if (is_vector<T_y>::value || is_vector<T_loc>::value) {
    using array_bool = Eigen::Array<bool, Eigen::Dynamic, 1>;
    cdf = forward_as<array_bool>(y_val < mu_val)
              .select(forward_as<T_partials_array>(exp_scaled_diff * 0.5),
                      1.0 - 0.5 / exp_scaled_diff)
              .prod();
    rep_deriv = forward_as<T_rep_deriv>(
        forward_as<array_bool>(y_val < mu_val)
            .select((cdf * inv_sigma),
                    forward_as<T_partials_array>(cdf * inv_sigma
                                                 / (2 * exp_scaled_diff - 1))));
  } else {
    if (is_vector<T_scale>::value) {
      cdf = forward_as<bool>(y_val < mu_val)
                ? forward_as<T_partials_array>(exp_scaled_diff * 0.5).prod()
                : forward_as<T_partials_array>(1.0 - 0.5 / exp_scaled_diff)
                      .prod();
    } else {
      cdf = forward_as<bool>(y_val < mu_val)
                ? forward_as<T_partials_return>(exp_scaled_diff * 0.5)
                : forward_as<T_partials_return>(1.0 - 0.5 / exp_scaled_diff);
    }
    if (forward_as<bool>(y_val < mu_val)) {
      rep_deriv = cdf * inv_sigma;
    } else {
      rep_deriv = cdf * inv_sigma / (2 * exp_scaled_diff - 1);
    }
  }

  if (!is_constant_all<T_y>::value) {
    partials<0>(ops_partials) = rep_deriv;
  }
  if (!is_constant_all<T_loc>::value) {
    partials<1>(ops_partials) = -rep_deriv;
  }
  if (!is_constant_all<T_scale>::value) {
    partials<2>(ops_partials) = -rep_deriv * scaled_diff;
  }
  return ops_partials.build(cdf);
}

}  // namespace math
}  // namespace stan
#endif
