#ifndef STAN_MATH_PRIM_PROB_DOUBLE_EXPONENTIAL_CDF_HPP
#define STAN_MATH_PRIM_PROB_DOUBLE_EXPONENTIAL_CDF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/exp.hpp>
#include <stan/math/prim/fun/inv.hpp>
#include <stan/math/prim/fun/max_size.hpp>
#include <stan/math/prim/fun/size_zero.hpp>
#include <stan/math/prim/fun/value_of.hpp>
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
template <typename T_y, typename T_loc, typename T_scale>
return_type_t<T_y, T_loc, T_scale> double_exponential_cdf(
    const T_y& y, const T_loc& mu, const T_scale& sigma) {
  using T_partials_return = partials_return_t<T_y, T_loc, T_scale>;
  using T_partials_array = Eigen::Array<T_partials_return, Eigen::Dynamic, 1>;
  using T_rep_deriv
      = std::conditional_t<is_vector<T_y>::value || is_vector<T_loc>::value
                               || is_vector<T_scale>::value,
                           T_partials_array, T_partials_return>;
  using std::exp;
  using T_y_ref = ref_type_t<T_y>;
  using T_mu_ref = ref_type_t<T_loc>;
  using T_sigma_ref = ref_type_t<T_scale>;
  static const char* function = "double_exponential_cdf";
  T_y_ref y_ref = y;
  T_mu_ref mu_ref = mu;
  T_sigma_ref sigma_ref = sigma;
  check_not_nan(function, "Random variable", y_ref);
  check_finite(function, "Location parameter", mu_ref);
  check_positive_finite(function, "Scale parameter", sigma_ref);

  if (size_zero(y, mu, sigma)) {
    return 1.0;
  }

  T_partials_return cdf(1.0);
  operands_and_partials<T_y_ref, T_mu_ref, T_sigma_ref> ops_partials(
      y_ref, mu_ref, sigma_ref);

  const auto& y_col = as_column_vector_or_scalar(y_ref);
  const auto& mu_col = as_column_vector_or_scalar(mu_ref);
  const auto& sigma_col = as_column_vector_or_scalar(sigma_ref);

  const auto& y_arr = as_array_or_scalar(y_col);
  const auto& mu_arr = as_array_or_scalar(mu_col);
  const auto& sigma_arr = as_array_or_scalar(sigma_col);

  ref_type_t<decltype(value_of(y_arr))> y_val = value_of(y_arr);
  ref_type_t<decltype(value_of(mu_arr))> mu_val = value_of(mu_arr);
  ref_type_t<decltype(value_of(sigma_arr))> sigma_val = value_of(sigma_arr);

  const auto& inv_sigma = to_ref(inv(sigma_val));
  const auto& scaled_diff = to_ref_if<!is_constant_all<T_scale>::value>(
      (y_val - mu_val) * inv_sigma);
  const auto& exp_scaled_diff = to_ref(exp(scaled_diff));

  size_t N = max_size(y, mu, sigma);
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
                    forward_as<T_partials_array>(cdf * inv_sigma / (2 * exp_scaled_diff - 1))));
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

  if (!is_constant_all<T_loc>::value) {
    if (is_vector<T_loc>::value) {
      ops_partials.edge2_.partials_ = -forward_as<T_partials_array>(rep_deriv);
    } else {
      ops_partials.edge2_.partials_[0] = -sum(rep_deriv);
    }
  }
  if (!is_constant_all<T_scale>::value) {
    if (is_vector<T_scale>::value) {
      ops_partials.edge3_.partials_
          = -forward_as<T_partials_array>(rep_deriv * scaled_diff);
    } else {
      ops_partials.edge3_.partials_[0] = -sum(rep_deriv) * sum(scaled_diff);
    }
  }
  if (!is_constant_all<T_y>::value) {
    if (is_vector<T_y>::value) {
      ops_partials.edge1_.partials_
          = std::move(forward_as<T_partials_array>(rep_deriv));
    } else {
      ops_partials.edge1_.partials_[0] = sum(rep_deriv);
    }
  }

  //  scalar_seq_view<T_y_ref> y_vec(y_ref);
  //  scalar_seq_view<T_mu_ref> mu_vec(mu_ref);
  //  scalar_seq_view<T_sigma_ref> sigma_vec(sigma_ref);
  //  size_t size_sigma = stan::math::size(sigma);
  //  size_t N = max_size(y, mu, sigma);

  //  VectorBuilder<true, T_partials_return, T_scale> inv_sigma(size_sigma);
  //  for (size_t i = 0; i < size_sigma; i++) {
  //    inv_sigma[i] = inv(value_of(sigma_vec[i]));
  //  }

  //  VectorBuilder<true, T_partials_return, T_y, T_loc, T_scale>
  //  scaled_diff(N); VectorBuilder<true, T_partials_return, T_y, T_loc,
  //  T_scale> exp_scaled_diff(
  //      N);

  //  for (size_t n = 0; n < N; n++) {
  //    const T_partials_return y_dbl = value_of(y_vec[n]);
  //    const T_partials_return mu_dbl = value_of(mu_vec[n]);
  //    scaled_diff[n] = (y_dbl - mu_dbl) * inv_sigma[n];
  //    exp_scaled_diff[n] = exp(scaled_diff[n]);

  //    if (y_dbl < mu_dbl) {
  //      cdf *= exp_scaled_diff[n] * 0.5;
  //    } else {
  //      cdf *= 1.0 - 0.5 / exp_scaled_diff[n];
  //    }
  //  }

  //  for (size_t n = 0; n < N; n++) {
  //    const T_partials_return y_dbl = value_of(y_vec[n]);
  //    const T_partials_return mu_dbl = value_of(mu_vec[n]);
  //    const T_partials_return sigma_dbl = value_of(sigma_vec[n]);

  //    const T_partials_return rep_deriv
  //        = y_dbl < mu_dbl ? cdf * inv_sigma[n]
  //                         : cdf * inv_sigma[n] / (2 * exp_scaled_diff[n] -
  //                         1);

  //    if (!is_constant_all<T_y>::value) {
  //      ops_partials.edge1_.partials_[n] += rep_deriv;
  //    }
  //    if (!is_constant_all<T_loc>::value) {
  //      ops_partials.edge2_.partials_[n] -= rep_deriv;
  //    }
  //    if (!is_constant_all<T_scale>::value) {
  //      ops_partials.edge3_.partials_[n] -= rep_deriv * scaled_diff[n];
  //    }
  //  }
  return ops_partials.build(cdf);
}

}  // namespace math
}  // namespace stan
#endif
