#ifndef STAN_MATH_PRIM_SCAL_PROB_CAUCHY_LCCDF_HPP
#define STAN_MATH_PRIM_SCAL_PROB_CAUCHY_LCCDF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/scal/err/check_consistent_sizes.hpp>
#include <stan/math/prim/scal/err/check_finite.hpp>
#include <stan/math/prim/scal/err/check_not_nan.hpp>
#include <stan/math/prim/scal/err/check_positive_finite.hpp>
#include <stan/math/prim/scal/fun/size_zero.hpp>
#include <stan/math/prim/scal/fun/constants.hpp>
#include <stan/math/prim/scal/fun/value_of.hpp>
#include <cmath>
#include <utility>

namespace stan {
namespace math {

/**
 * Returns the cauchy log complementary cumulative distribution function
 * for the given location, and scale. If given containers of matching sizes
 * returns the log sum of probabilities.
 *
 * @tparam T_y type of real parameter
 * @tparam T_loc type of location parameter
 * @tparam T_scale type of scale parameter
 * @param y real parameter
 * @param mu location parameter
 * @param sigma scale parameter
 * @return log probability or log sum of probabilities
 * @throw std::domain_error if sigma is nonpositive or y, mu are nan
 * @throw std::invalid_argument if container sizes mismatch
 */
template <typename T_y, typename T_loc, typename T_scale>
inline auto cauchy_lccdf(T_y&& y, T_loc&& mu, T_scale&& sigma) {
  using T_partials = partials_return_t<T_y, T_loc, T_scale>;
  T_partials ccdf_log(0.0);
  using std::atan;
  using std::log;

  static const char* function = "cauchy_lccdf";
  check_not_nan(function, "Random variable", y);
  check_finite(function, "Location parameter", mu);
  check_positive_finite(function, "Scale parameter", sigma);
  check_consistent_sizes(function, "Random variable", y, "Location parameter",
                         mu, "Scale Parameter", sigma);

  const scalar_seq_view<T_y> y_vec(y);
  const scalar_seq_view<T_loc> mu_vec(mu);
  const scalar_seq_view<T_scale> sigma_vec(sigma);
  const auto N = max_size(y, mu, sigma);
  operands_and_partials<T_y, T_loc, T_scale> ops_partials(y, mu, sigma);
  if (size_zero(y, mu, sigma)) {
    return ops_partials.build(ccdf_log);
  }

  for (size_t n = 0; n < N; n++) {
    const auto y_dbl = value_of(y_vec[n]);
    const auto mu_dbl = value_of(mu_vec[n]);
    const auto sigma_inv_dbl = 1.0 / value_of(sigma_vec[n]);
    const auto sigma_dbl = value_of(sigma_vec[n]);
    const auto z = (y_dbl - mu_dbl) * sigma_inv_dbl;

    const auto Pn = 0.5 - atan(z) / pi();
    ccdf_log += log(Pn);

    const auto rep_deriv = 1.0 / (Pn * pi() * (z * z * sigma_dbl + sigma_dbl));
    if (!is_constant_all<T_y>::value) {
      ops_partials.edge1_.partials_[n] -= rep_deriv;
    }
    if (!is_constant_all<T_loc>::value) {
      ops_partials.edge2_.partials_[n] += rep_deriv;
    }
    if (!is_constant_all<T_scale>::value) {
      ops_partials.edge3_.partials_[n] += rep_deriv * z;
    }
  }
  return ops_partials.build(ccdf_log);
}

}  // namespace math
}  // namespace stan
#endif
