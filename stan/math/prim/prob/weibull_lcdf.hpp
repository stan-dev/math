#ifndef STAN_MATH_PRIM_PROB_WEIBULL_LCDF_HPP
#define STAN_MATH_PRIM_PROB_WEIBULL_LCDF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/exp.hpp>
#include <stan/math/prim/fun/log.hpp>
#include <stan/math/prim/fun/size_zero.hpp>
#include <stan/math/prim/fun/max_size.hpp>
#include <stan/math/prim/fun/value_of.hpp>
#include <stan/math/prim/functor/operands_and_partials.hpp>
#include <cmath>

namespace stan {
namespace math {

/** \ingroup prob_dists
 * Returns the Weibull log cumulative distribution function for the given
 * location and scale. Given containers of matching sizes, returns the
 * log sum of probabilities.
 *
 * @tparam T_y type of real parameter
 * @tparam T_shape type of shape parameter
 * @tparam T_scale type of scale parameter
 * @param y real parameter
 * @param alpha shape parameter
 * @param sigma scale parameter
 * @return log probability or log sum of probabilities
 * @throw std::domain_error if y is negative, alpha sigma is nonpositive
 */
template <typename T_y, typename T_shape, typename T_scale>
return_type_t<T_y, T_shape, T_scale> weibull_lcdf(const T_y& y,
                                                  const T_shape& alpha,
                                                  const T_scale& sigma) {
  using T_partials_return = partials_return_t<T_y, T_shape, T_scale>;
  using std::exp;
  using std::log;
  using std::pow;
  static const char* function = "weibull_lcdf";
  check_nonnegative(function, "Random variable", y);
  check_positive_finite(function, "Shape parameter", alpha);
  check_positive_finite(function, "Scale parameter", sigma);

  if (size_zero(y, alpha, sigma)) {
    return 0.0;
  }

  T_partials_return cdf_log(0.0);
  operands_and_partials<T_y, T_shape, T_scale> ops_partials(y, alpha, sigma);

  scalar_seq_view<T_y> y_vec(y);
  scalar_seq_view<T_scale> sigma_vec(sigma);
  scalar_seq_view<T_shape> alpha_vec(alpha);
  size_t N = max_size(y, sigma, alpha);
  for (size_t n = 0; n < N; n++) {
    const T_partials_return y_dbl = value_of(y_vec[n]);
    const T_partials_return sigma_dbl = value_of(sigma_vec[n]);
    const T_partials_return alpha_dbl = value_of(alpha_vec[n]);
    const T_partials_return pow_n = pow(y_dbl / sigma_dbl, alpha_dbl);
    const T_partials_return exp_n = exp(-pow_n);
    const T_partials_return cdf_n = 1.0 - exp_n;

    cdf_log += log(cdf_n);

    const T_partials_return rep_deriv = pow_n / (1.0 / exp_n - 1.0);
    if (!is_constant_all<T_y>::value) {
      ops_partials.edge1_.partials_[n] += rep_deriv * alpha_dbl / y_dbl;
    }
    if (!is_constant_all<T_shape>::value) {
      ops_partials.edge2_.partials_[n] += rep_deriv * log(y_dbl / sigma_dbl);
    }
    if (!is_constant_all<T_scale>::value) {
      ops_partials.edge3_.partials_[n] -= rep_deriv * alpha_dbl / sigma_dbl;
    }
  }
  return ops_partials.build(cdf_log);
}

}  // namespace math
}  // namespace stan
#endif
