#ifndef STAN_MATH_PRIM_PROB_EXPONENTIAL_CDF_HPP
#define STAN_MATH_PRIM_PROB_EXPONENTIAL_CDF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/size_zero.hpp>
#include <stan/math/prim/fun/value_of.hpp>
#include <cmath>

namespace stan {
namespace math {

/** \ingroup prob_dists
 * Calculates the exponential cumulative distribution function for
 * the given y and beta.
 *
 * Inverse scale parameter must be greater than 0.
 * y must be greater than or equal to 0.
 *
 * @param y A scalar variable.
 * @param beta Inverse scale parameter.
 * @tparam T_y Type of scalar.
 * @tparam T_inv_scale Type of inverse scale.
 */
template <typename T_y, typename T_inv_scale>
return_type_t<T_y, T_inv_scale> exponential_cdf(const T_y& y,
                                                const T_inv_scale& beta) {
  using T_partials_return = partials_return_t<T_y, T_inv_scale>;

  static const char* function = "exponential_cdf";

  using std::exp;

  T_partials_return cdf(1.0);
  if (size_zero(y, beta)) {
    return cdf;
  }

  check_not_nan(function, "Random variable", y);
  check_nonnegative(function, "Random variable", y);
  check_positive_finite(function, "Inverse scale parameter", beta);

  operands_and_partials<T_y, T_inv_scale> ops_partials(y, beta);

  scalar_seq_view<T_y> y_vec(y);
  scalar_seq_view<T_inv_scale> beta_vec(beta);
  size_t N = max_size(y, beta);
  for (size_t n = 0; n < N; n++) {
    const T_partials_return beta_dbl = value_of(beta_vec[n]);
    const T_partials_return y_dbl = value_of(y_vec[n]);
    const T_partials_return one_m_exp = 1.0 - exp(-beta_dbl * y_dbl);

    cdf *= one_m_exp;
  }

  for (size_t n = 0; n < N; n++) {
    const T_partials_return beta_dbl = value_of(beta_vec[n]);
    const T_partials_return y_dbl = value_of(y_vec[n]);
    const T_partials_return one_m_exp = 1.0 - exp(-beta_dbl * y_dbl);

    T_partials_return rep_deriv = exp(-beta_dbl * y_dbl) / one_m_exp;
    if (!is_constant_all<T_y>::value) {
      ops_partials.edge1_.partials_[n] += rep_deriv * beta_dbl * cdf;
    }
    if (!is_constant_all<T_inv_scale>::value) {
      ops_partials.edge2_.partials_[n] += rep_deriv * y_dbl * cdf;
    }
  }
  return ops_partials.build(cdf);
}

}  // namespace math
}  // namespace stan
#endif
