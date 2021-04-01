#ifndef STAN_MATH_PRIM_PROB_DOUBLE_EXPONENTIAL_LCDF_HPP
#define STAN_MATH_PRIM_PROB_DOUBLE_EXPONENTIAL_LCDF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/exp.hpp>
#include <stan/math/prim/fun/inv.hpp>
#include <stan/math/prim/fun/log1m.hpp>
#include <stan/math/prim/fun/max_size.hpp>
#include <stan/math/prim/fun/scalar_seq_view.hpp>
#include <stan/math/prim/fun/size_zero.hpp>
#include <stan/math/prim/fun/value_of.hpp>
#include <stan/math/prim/functor/operands_and_partials.hpp>
#include <cmath>

namespace stan {
namespace math {

/** \ingroup prob_dists
 * Returns the double exponential log cumulative density function. Given
 * containers of matching sizes, returns the log sum of probabilities.
 *
 * @tparam T_y type of real parameter.
 * @tparam T_loc type of location parameter.
 * @tparam T_scale type of scale parameter.
 * @param y real parameter
 * @param mu location parameter
 * @param sigma scale parameter
 * @return log probability or log sum of probabilities
 * @throw std::domain_error if y is nan, mu is infinite, or sigma is nonpositive
 * @throw std::invalid_argument if container sizes mismatch
 */
template <typename T_y, typename T_loc, typename T_scale,
          require_all_not_nonscalar_prim_or_rev_kernel_expression_t<
              T_y, T_loc, T_scale>* = nullptr>
return_type_t<T_y, T_loc, T_scale> double_exponential_lcdf(
    const T_y& y, const T_loc& mu, const T_scale& sigma) {
  using T_partials_return = partials_return_t<T_y, T_loc, T_scale>;
  using std::exp;
  using std::log;
  using T_y_ref = ref_type_t<T_y>;
  using T_mu_ref = ref_type_t<T_loc>;
  using T_sigma_ref = ref_type_t<T_scale>;
  static const char* function = "double_exponential_lcdf";
  check_consistent_sizes(function, "Random variable", y, "Location parameter",
                         mu, "Scale Parameter", sigma);
  T_y_ref y_ref = y;
  T_mu_ref mu_ref = mu;
  T_sigma_ref sigma_ref = sigma;
  check_not_nan(function, "Random variable", y_ref);
  check_finite(function, "Location parameter", mu_ref);
  check_positive_finite(function, "Scale parameter", sigma_ref);

  if (size_zero(y, mu, sigma)) {
    return 0;
  }

  T_partials_return cdf_log(0.0);
  operands_and_partials<T_y_ref, T_mu_ref, T_sigma_ref> ops_partials(
      y_ref, mu_ref, sigma_ref);

  scalar_seq_view<T_y_ref> y_vec(y_ref);
  scalar_seq_view<T_mu_ref> mu_vec(mu_ref);
  scalar_seq_view<T_sigma_ref> sigma_vec(sigma_ref);
  size_t size_sigma = stan::math::size(sigma);
  size_t N = max_size(y, mu, sigma);

  VectorBuilder<true, T_partials_return, T_scale> inv_sigma(size_sigma);
  for (size_t i = 0; i < size_sigma; i++) {
    inv_sigma[i] = inv(sigma_vec.val(i));
  }

  for (size_t n = 0; n < N; n++) {
    const T_partials_return y_dbl = y_vec.val(n);
    const T_partials_return mu_dbl = mu_vec.val(n);
    const T_partials_return scaled_diff = (y_dbl - mu_dbl) * inv_sigma[n];

    const T_partials_return rep_deriv
        = y_dbl < mu_dbl ? inv_sigma[n]
                         : inv_sigma[n] * inv(2 * exp(scaled_diff) - 1);

    if (y_dbl < mu_dbl) {
      cdf_log += LOG_HALF + scaled_diff;
    } else {
      cdf_log += log1m(0.5 * exp(-scaled_diff));
    }

    if (!is_constant_all<T_y>::value) {
      ops_partials.edge1_.partials_[n] += rep_deriv;
    }
    if (!is_constant_all<T_loc>::value) {
      ops_partials.edge2_.partials_[n] -= rep_deriv;
    }
    if (!is_constant_all<T_scale>::value) {
      ops_partials.edge3_.partials_[n] -= rep_deriv * scaled_diff;
    }
  }
  return ops_partials.build(cdf_log);
}
}  // namespace math
}  // namespace stan
#endif
