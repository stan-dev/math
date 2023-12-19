#ifndef STAN_MATH_PRIM_PROB_NORMAL_SUFFICIENT_LPDF_HPP
#define STAN_MATH_PRIM_PROB_NORMAL_SUFFICIENT_LPDF_HPP

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
#include <stan/math/prim/prob/normal_lpdf.hpp>
#include <stan/math/prim/functor/partials_propagator.hpp>
#include <cmath>

namespace stan {
namespace math {

/** \ingroup prob_dists
 * The log of the normal density for the specified scalar(s) given
 * the specified mean(s) and deviation(s).
 * y, s_squared, mu, or sigma can each be either
 * a scalar, a std vector or Eigen vector.
 * n can be either a single int or an std vector of ints.
 * Any vector inputs must be the same length.
 *
 * <p>The result log probability is defined to be the sum of the
 * log probabilities for each observation/mean/deviation triple.
 *
 * @tparam T_y type of sample average parameter
 * @tparam T_s type of sample squared errors parameter
 * @tparam T_n type of sample size parameter
 * @tparam T_loc type of location parameter
 * @tparam T_scale type of scale parameter
 *
 * @param y_bar (Sequence of) scalar(s) (sample average(s)).
 * @param s_squared (Sequence of) sum(s) of sample squared errors
 * @param n_obs (Sequence of) sample size(s)
 * @param mu (Sequence of) location parameter(s)
 * for the normal distribution.
 * @param sigma (Sequence of) scale parameters for the normal
 * distribution.
 * @return The log of the product of the densities.
 * @throw std::domain_error if either n or sigma are not positive,
 * if s_squared is negative or if any parameter is not finite.
 */
template <bool propto, typename T_y, typename T_s, typename T_n, typename T_loc,
          typename T_scale>
return_type_t<T_y, T_s, T_loc, T_scale> normal_sufficient_lpdf(
    const T_y& y_bar, const T_s& s_squared, const T_n& n_obs, const T_loc& mu,
    const T_scale& sigma) {
  using T_partials_return = partials_return_t<T_y, T_s, T_n, T_loc, T_scale>;
  using T_y_ref = ref_type_if_not_constant_t<T_y>;
  using T_s_ref = ref_type_if_not_constant_t<T_s>;
  using T_n_ref = ref_type_if_not_constant_t<T_n>;
  using T_mu_ref = ref_type_if_not_constant_t<T_loc>;
  using T_sigma_ref = ref_type_if_not_constant_t<T_scale>;
  static const char* function = "normal_sufficient_lpdf";
  check_consistent_sizes(function, "Location parameter sufficient statistic",
                         y_bar, "Scale parameter sufficient statistic",
                         s_squared, "Number of observations", n_obs,
                         "Location parameter", mu, "Scale parameter", sigma);

  T_y_ref y_ref = y_bar;
  T_s_ref s_squared_ref = s_squared;
  T_n_ref n_obs_ref = n_obs;
  T_mu_ref mu_ref = mu;
  T_sigma_ref sigma_ref = sigma;

  decltype(auto) y_val = to_ref(as_value_column_array_or_scalar(y_ref));
  decltype(auto) s_squared_val
      = to_ref(as_value_column_array_or_scalar(s_squared_ref));
  decltype(auto) n_obs_val_int
      = to_ref(as_value_column_array_or_scalar(n_obs_ref));
  decltype(auto) n_obs_val = to_ref(
      promote_scalar<double>(as_value_column_array_or_scalar(n_obs_ref)));
  decltype(auto) mu_val = to_ref(as_value_column_array_or_scalar(mu_ref));
  decltype(auto) sigma_val = to_ref(as_value_column_array_or_scalar(sigma_ref));

  check_finite(function, "Location parameter sufficient statistic", y_val);
  check_finite(function, "Scale parameter sufficient statistic", s_squared_val);
  check_nonnegative(function, "Scale parameter sufficient statistic",
                    s_squared_val);
  check_positive_finite(function, "Number of observations", n_obs_val);
  check_finite(function, "Location parameter", mu_val);
  check_positive_finite(function, "Scale parameter", sigma_val);

  if (size_zero(y_bar, s_squared, n_obs, mu, sigma)) {
    return 0.0;
  }
  if (!include_summand<propto, T_y, T_s, T_loc, T_scale>::value) {
    return 0.0;
  }

  const auto& sigma_squared
      = to_ref_if<!is_constant_all<T_y, T_loc, T_s, T_scale>::value>(
          square(sigma_val));
  const auto& diff = to_ref(mu_val - y_val);
  const auto& cons_expr = to_ref_if<!is_constant_all<T_scale>::value>(
      s_squared_val + n_obs_val * diff * diff);

  size_t N = max_size(y_bar, s_squared, n_obs, mu, sigma);
  T_partials_return logp = -sum(cons_expr / (2 * sigma_squared));
  if (include_summand<propto>::value) {
    logp += NEG_LOG_SQRT_TWO_PI * sum(n_obs_val) * N / math::size(n_obs);
  }
  if (include_summand<propto, T_scale>::value) {
    logp -= sum(n_obs_val * log(sigma_val)) * N / max_size(n_obs, sigma);
  }

  auto ops_partials
      = make_partials_propagator(y_ref, s_squared_ref, mu_ref, sigma_ref);
  if (!is_constant_all<T_y, T_loc>::value) {
    auto common_derivative = to_ref_if<(!is_constant_all<T_loc>::value
                                        && !is_constant_all<T_y>::value)>(
        N / max_size(y_bar, mu, n_obs, sigma) * n_obs_val / sigma_squared
        * diff);
    if (!is_constant_all<T_loc>::value) {
      partials<2>(ops_partials) = -common_derivative;
    }
    if (!is_constant_all<T_y>::value) {
      partials<0>(ops_partials) = std::move(common_derivative);
    }
  }
  if (!is_constant_all<T_s>::value) {
    using T_sigma_value_scalar = scalar_type_t<decltype(sigma_val)>;
    using T_sigma_value_vector
        = Eigen::Array<T_sigma_value_scalar, Eigen::Dynamic, 1>;
    if (is_vector<T_scale>::value) {
      edge<1>(ops_partials).partials_
          = -0.5 / forward_as<T_sigma_value_vector>(sigma_squared);
    } else {
      if (is_vector<T_s>::value) {
        partials<1>(ops_partials) = T_sigma_value_vector::Constant(
            N, -0.5 / forward_as<T_sigma_value_scalar>(sigma_squared));
      } else {
        forward_as<internal::broadcast_array<T_partials_return>>(
            partials<1>(ops_partials))
            = -0.5 / sigma_squared * N / math::size(sigma);
      }
    }
  }
  if (!is_constant_all<T_scale>::value) {
    edge<3>(ops_partials).partials_
        = (cons_expr / sigma_squared - n_obs_val) / sigma_val;
  }
  return ops_partials.build(logp);
}

template <typename T_y, typename T_s, typename T_n, typename T_loc,
          typename T_scale>
inline return_type_t<T_y, T_s, T_loc, T_scale> normal_sufficient_lpdf(
    const T_y& y_bar, const T_s& s_squared, const T_n& n_obs, const T_loc& mu,
    const T_scale& sigma) {
  return normal_sufficient_lpdf<false>(y_bar, s_squared, n_obs, mu, sigma);
}

}  // namespace math
}  // namespace stan
#endif
