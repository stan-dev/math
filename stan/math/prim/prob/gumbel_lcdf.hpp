#ifndef STAN_MATH_PRIM_PROB_GUMBEL_LCDF_HPP
#define STAN_MATH_PRIM_PROB_GUMBEL_LCDF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/as_column_vector_or_scalar.hpp>
#include <stan/math/prim/fun/as_array_or_scalar.hpp>
#include <stan/math/prim/fun/as_value_column_array_or_scalar.hpp>
#include <stan/math/prim/fun/exp.hpp>
#include <stan/math/prim/fun/max_size.hpp>
#include <stan/math/prim/fun/size_zero.hpp>
#include <stan/math/prim/fun/to_ref.hpp>
#include <stan/math/prim/fun/value_of.hpp>
#include <stan/math/prim/functor/partials_propagator.hpp>
#include <cmath>

namespace stan {
namespace math {

/** \ingroup prob_dists
 * Returns the Gumbel log cumulative distribution for the given
 * location and scale. Given containers of matching sizes, returns the
 * log sum of probabilities.
 *
 * @tparam T_y type of real parameter
 * @tparam T_loc type of location parameter
 * @tparam T_scale type of scale parameter
 * @param y real parameter
 * @param mu location parameter
 * @param beta scale parameter
 * @return log probability or log sum of probabilities
 * @throw std::domain_error if y is nan, mu is infinite, or beta is nonpositive
 * @throw std::invalid_argument if container sizes mismatch
 */
template <typename T_y, typename T_loc, typename T_scale,
          require_all_not_nonscalar_prim_or_rev_kernel_expression_t<
              T_y, T_loc, T_scale>* = nullptr>
return_type_t<T_y, T_loc, T_scale> gumbel_lcdf(const T_y& y, const T_loc& mu,
                                               const T_scale& beta) {
  using T_partials_return = partials_return_t<T_y, T_loc, T_scale>;
  using T_y_ref = ref_type_if_t<!is_constant<T_y>::value, T_y>;
  using T_mu_ref = ref_type_if_t<!is_constant<T_loc>::value, T_loc>;
  using T_beta_ref = ref_type_if_t<!is_constant<T_scale>::value, T_scale>;
  static const char* function = "gumbel_lcdf";
  check_consistent_sizes(function, "Random variable", y, "Location parameter",
                         mu, "Scale parameter", beta);
  T_y_ref y_ref = y;
  T_mu_ref mu_ref = mu;
  T_beta_ref beta_ref = beta;

  decltype(auto) y_val = to_ref(as_value_column_array_or_scalar(y_ref));
  decltype(auto) mu_val = to_ref(as_value_column_array_or_scalar(mu_ref));
  decltype(auto) beta_val = to_ref(as_value_column_array_or_scalar(beta_ref));

  check_not_nan(function, "Random variable", y_val);
  check_finite(function, "Location parameter", mu_val);
  check_positive(function, "Scale parameter", beta_val);

  if (size_zero(y, mu, beta)) {
    return 0;
  }

  auto ops_partials = make_partials_propagator(y_ref, mu_ref, beta_ref);

  const auto& scaled_diff
      = to_ref_if<!is_constant_all<T_y, T_loc, T_scale>::value>((mu_val - y_val)
                                                                / beta_val);
  const auto& exp_scaled_diff
      = to_ref_if<!is_constant_all<T_y, T_loc, T_scale>::value>(
          exp(scaled_diff));
  T_partials_return cdf_log = -sum(exp_scaled_diff);

  if (!is_constant_all<T_y, T_loc, T_scale>::value) {
    const auto& rep_deriv = to_ref_if<!is_constant_all<T_loc>::value
                                          + !is_constant_all<T_scale>::value
                                          + !is_constant_all<T_y>::value
                                      >= 2>(exp_scaled_diff / beta_val);
    if (!is_constant_all<T_y>::value) {
      partials<0>(ops_partials) = rep_deriv;
    }
    if (!is_constant_all<T_loc>::value) {
      partials<1>(ops_partials) = -rep_deriv;
    }
    if (!is_constant_all<T_scale>::value) {
      partials<2>(ops_partials) = rep_deriv * scaled_diff;
    }
  }
  return ops_partials.build(cdf_log);
}

}  // namespace math
}  // namespace stan
#endif
