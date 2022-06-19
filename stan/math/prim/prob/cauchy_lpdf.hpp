#ifndef STAN_MATH_PRIM_PROB_CAUCHY_LPDF_HPP
#define STAN_MATH_PRIM_PROB_CAUCHY_LPDF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/as_column_vector_or_scalar.hpp>
#include <stan/math/prim/fun/as_array_or_scalar.hpp>
#include <stan/math/prim/fun/as_value_column_array_or_scalar.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/log.hpp>
#include <stan/math/prim/fun/log1p.hpp>
#include <stan/math/prim/fun/max_size.hpp>
#include <stan/math/prim/fun/size.hpp>
#include <stan/math/prim/fun/size_zero.hpp>
#include <stan/math/prim/fun/to_ref.hpp>
#include <stan/math/prim/fun/value_of.hpp>
#include <stan/math/prim/functor/operands_and_partials.hpp>
#include <cmath>

namespace stan {
namespace math {

/** \ingroup prob_dists
 * The log of the Cauchy density for the specified scalar(s) given
 * the specified location parameter(s) and scale parameter(s). y,
 * mu, or sigma can each either be scalar a vector.  Any vector
 * inputs must be the same length.
 *
 * <p> The result log probability is defined to be the sum of
 * the log probabilities for each observation/mu/sigma triple.
 *
 * @tparam T_y type of scalar outcome
 * @tparam T_loc type of location
 * @tparam T_scale type of scale
 * @param y (Sequence of) scalar(s).
 * @param mu (Sequence of) location(s).
 * @param sigma (Sequence of) scale(s).
 * @return The log of the product of densities.
 */
template <bool propto, typename T_y, typename T_loc, typename T_scale,
          require_all_not_nonscalar_prim_or_rev_kernel_expression_t<
              T_y, T_loc, T_scale>* = nullptr>
return_type_t<T_y, T_loc, T_scale> cauchy_lpdf(const T_y& y, const T_loc& mu,
                                               const T_scale& sigma) {
  using T_partials_return = partials_return_t<T_y, T_loc, T_scale>;
  using std::log;
  using T_y_ref = ref_type_if_t<!is_constant<T_y>::value, T_y>;
  using T_mu_ref = ref_type_if_t<!is_constant<T_loc>::value, T_loc>;
  using T_sigma_ref = ref_type_if_t<!is_constant<T_scale>::value, T_scale>;
  static const char* function = "cauchy_lpdf";
  check_consistent_sizes(function, "Random variable", y, "Location parameter",
                         mu, "Scale parameter", sigma);
  T_y_ref y_ref = y;
  T_mu_ref mu_ref = mu;
  T_sigma_ref sigma_ref = sigma;

  if (size_zero(y, mu, sigma)) {
    return 0.0;
  }
  if (!include_summand<propto, T_y, T_loc, T_scale>::value) {
    return 0.0;
  }

  T_partials_return logp(0.0);
  operands_and_partials<T_y_ref, T_mu_ref, T_sigma_ref> ops_partials(
      y_ref, mu_ref, sigma_ref);

  decltype(auto) y_val = to_ref(as_value_column_array_or_scalar(y_ref));
  decltype(auto) mu_val = to_ref(as_value_column_array_or_scalar(mu_ref));
  decltype(auto) sigma_val = to_ref(as_value_column_array_or_scalar(sigma_ref));
  check_not_nan(function, "Random variable", y_val);
  check_finite(function, "Location parameter", mu_val);
  check_positive_finite(function, "Scale parameter", sigma_val);

  size_t N = max_size(y, mu, sigma);

  const auto& inv_sigma
      = to_ref_if<!is_constant_all<T_scale>::value>(inv(sigma_val));
  const auto& y_minus_mu
      = to_ref_if<!is_constant_all<T_y, T_loc, T_scale>::value>(y_val - mu_val);

  logp -= sum(log1p(square(y_minus_mu * inv_sigma)));
  if (include_summand<propto>::value) {
    logp -= N * LOG_PI;
  }
  if (include_summand<propto, T_scale>::value) {
    logp -= sum(log(sigma_val)) * N / math::size(sigma);
  }

  if (!is_constant_all<T_y, T_loc, T_scale>::value) {
    const auto& sigma_squared
        = to_ref_if<!is_constant_all<T_scale>::value>(square(sigma_val));
    const auto& y_minus_mu_squared
        = to_ref_if<!is_constant_all<T_scale>::value>(square(y_minus_mu));
    if (!is_constant_all<T_y, T_loc>::value) {
      auto mu_deriv = to_ref_if<(!is_constant_all<T_y>::value
                                 && !is_constant_all<T_loc>::value)>(
          2 * y_minus_mu / (sigma_squared + y_minus_mu_squared));
      if (!is_constant_all<T_y>::value) {
        if (is_vector<T_y>::value) {
          ops_partials.edge1_.partials_ = -mu_deriv;
        } else {
          ops_partials.edge1_.partials_[0] = -sum(mu_deriv);
        }
      }
      if (!is_constant_all<T_loc>::value) {
        ops_partials.edge2_.partials_ = std::move(mu_deriv);
      }
    }
    if (!is_constant_all<T_scale>::value) {
      ops_partials.edge3_.partials_ = (y_minus_mu_squared - sigma_squared)
                                      * inv_sigma
                                      / (sigma_squared + y_minus_mu_squared);
    }
  }
  return ops_partials.build(logp);
}

template <typename T_y, typename T_loc, typename T_scale>
inline return_type_t<T_y, T_loc, T_scale> cauchy_lpdf(const T_y& y,
                                                      const T_loc& mu,
                                                      const T_scale& sigma) {
  return cauchy_lpdf<false>(y, mu, sigma);
}

}  // namespace math
}  // namespace stan
#endif
