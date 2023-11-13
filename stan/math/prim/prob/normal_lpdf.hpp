#ifndef STAN_MATH_PRIM_PROB_NORMAL_LPDF_HPP
#define STAN_MATH_PRIM_PROB_NORMAL_LPDF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/as_column_vector_or_scalar.hpp>
#include <stan/math/prim/fun/as_array_or_scalar.hpp>
#include <stan/math/prim/fun/as_value_column_array_or_scalar.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/log.hpp>
#include <stan/math/prim/fun/max_size.hpp>
#include <stan/math/prim/fun/size.hpp>
#include <stan/math/prim/fun/size_zero.hpp>
#include <stan/math/prim/fun/to_ref.hpp>
#include <stan/math/prim/fun/value_of.hpp>
#include <stan/math/prim/functor/partials_propagator.hpp>
#include <cmath>

namespace stan {
namespace math {

/** \ingroup prob_dists
 * The log of the normal density for the specified scalar(s) given
 * the specified mean(s) and deviation(s). y, mu, or sigma can
 * each be either a scalar or a vector. Any vector inputs
 * must be the same length.
 *
 * <p>The result log probability is defined to be the sum of the
 * log probabilities for each observation/mean/deviation triple.
 *
 * @tparam T_y type of scalar
 * @tparam T_loc type of location parameter
 * @tparam T_scale type of scale parameter
 * @param y (Sequence of) scalar(s).
 * @param mu (Sequence of) location parameter(s)
 * for the normal distribution.
 * @param sigma (Sequence of) scale parameters for the normal distribution.
 * @return The log of the product of the densities.
 * @throw std::domain_error if the scale is not positive.
 */
template <bool propto, typename T_y, typename T_loc, typename T_scale,
          require_all_not_nonscalar_prim_or_rev_kernel_expression_t<
              T_y, T_loc, T_scale>* = nullptr>
inline return_type_t<T_y, T_loc, T_scale> normal_lpdf(const T_y& y,
                                                      const T_loc& mu,
                                                      const T_scale& sigma) {
  using T_partials_return = partials_return_t<T_y, T_loc, T_scale>;
  using T_y_ref = ref_type_if_not_constant_t<T_y>;
  using T_mu_ref = ref_type_if_not_constant_t<T_loc>;
  using T_sigma_ref = ref_type_if_not_constant_t<T_scale>;
  static const char* function = "normal_lpdf";
  check_consistent_sizes(function, "Random variable", y, "Location parameter",
                         mu, "Scale parameter", sigma);
  T_y_ref y_ref = y;
  T_mu_ref mu_ref = mu;
  T_sigma_ref sigma_ref = sigma;

  decltype(auto) y_val = to_ref(as_value_column_array_or_scalar(y_ref));
  decltype(auto) mu_val = to_ref(as_value_column_array_or_scalar(mu_ref));
  decltype(auto) sigma_val = to_ref(as_value_column_array_or_scalar(sigma_ref));

  check_not_nan(function, "Random variable", y_val);
  check_finite(function, "Location parameter", mu_val);
  check_positive(function, "Scale parameter", sigma_val);

  if (size_zero(y, mu, sigma)) {
    return 0.0;
  }
  if (!include_summand<propto, T_y, T_loc, T_scale>::value) {
    return 0.0;
  }

  auto ops_partials = make_partials_propagator(y_ref, mu_ref, sigma_ref);

  const auto& inv_sigma
      = to_ref_if<!is_constant_all<T_y, T_scale, T_loc>::value>(inv(sigma_val));
  const auto& y_scaled = to_ref((y_val - mu_val) * inv_sigma);
  const auto& y_scaled_sq
      = to_ref_if<!is_constant_all<T_scale>::value>(y_scaled * y_scaled);

  size_t N = max_size(y, mu, sigma);
  T_partials_return logp = -0.5 * sum(y_scaled_sq);
  if (include_summand<propto>::value) {
    logp += NEG_LOG_SQRT_TWO_PI * N;
  }
  if (include_summand<propto, T_scale>::value) {
    logp -= sum(log(sigma_val)) * N / math::size(sigma);
  }

  if (!is_constant_all<T_y, T_scale, T_loc>::value) {
    auto scaled_diff = to_ref_if<!is_constant_all<T_y>::value
                                     + !is_constant_all<T_scale>::value
                                     + !is_constant_all<T_loc>::value
                                 >= 2>(inv_sigma * y_scaled);
    if (!is_constant_all<T_y>::value) {
      partials<0>(ops_partials) = -scaled_diff;
    }
    if (!is_constant_all<T_scale>::value) {
      partials<2>(ops_partials) = inv_sigma * y_scaled_sq - inv_sigma;
    }
    if (!is_constant_all<T_loc>::value) {
      partials<1>(ops_partials) = std::move(scaled_diff);
    }
  }
  return ops_partials.build(logp);
}

template <typename T_y, typename T_loc, typename T_scale>
inline return_type_t<T_y, T_loc, T_scale> normal_lpdf(const T_y& y,
                                                      const T_loc& mu,
                                                      const T_scale& sigma) {
  return normal_lpdf<false>(y, mu, sigma);
}

}  // namespace math
}  // namespace stan
#endif
