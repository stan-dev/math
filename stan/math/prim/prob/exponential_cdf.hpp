#ifndef STAN_MATH_PRIM_PROB_EXPONENTIAL_CDF_HPP
#define STAN_MATH_PRIM_PROB_EXPONENTIAL_CDF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/exp.hpp>
#include <stan/math/prim/fun/max_size.hpp>
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
 * @tparam T_y type of scalar
 * @tparam T_inv_scale type of inverse scale
 * @param y A scalar variable.
 * @param beta Inverse scale parameter.
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

  VectorBuilder<!is_constant_all<T_y, T_inv_scale>::value, T_partials_return,
                T_y, T_inv_scale>
      rep_deriv(N);
  for (size_t n = 0; n < N; n++) {
    const T_partials_return beta_dbl = value_of(beta_vec[n]);
    const T_partials_return y_dbl = value_of(y_vec[n]);
    const T_partials_return exp_val = exp(-beta_dbl * y_dbl);
    const T_partials_return one_m_exp = 1 - exp_val;

    if (!is_constant_all<T_y, T_inv_scale>::value) {
      rep_deriv[n] = exp_val / one_m_exp;
    }
    cdf *= one_m_exp;
  }

  for (size_t n = 0; n < N; n++) {
    if (!is_constant_all<T_y>::value) {
      ops_partials.edge1_.partials_[n]
          += value_of(beta_vec[n]) * rep_deriv[n] * cdf;
    }
    if (!is_constant_all<T_inv_scale>::value) {
      ops_partials.edge2_.partials_[n]
          += value_of(y_vec[n]) * rep_deriv[n] * cdf;
    }
  }
  return ops_partials.build(cdf);
}

}  // namespace math
}  // namespace stan
#endif
