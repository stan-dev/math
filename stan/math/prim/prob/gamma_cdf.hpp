#ifndef STAN_MATH_PRIM_PROB_GAMMA_CDF_HPP
#define STAN_MATH_PRIM_PROB_GAMMA_CDF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/digamma.hpp>
#include <stan/math/prim/fun/exp.hpp>
#include <stan/math/prim/fun/gamma_p.hpp>
#include <stan/math/prim/fun/grad_reg_inc_gamma.hpp>
#include <stan/math/prim/fun/max_size.hpp>
#include <stan/math/prim/fun/multiply_log.hpp>
#include <stan/math/prim/fun/size.hpp>
#include <stan/math/prim/fun/size_zero.hpp>
#include <stan/math/prim/fun/tgamma.hpp>
#include <stan/math/prim/fun/value_of.hpp>
#include <boost/random/gamma_distribution.hpp>
#include <boost/random/variate_generator.hpp>
#include <cmath>

namespace stan {
namespace math {

/** \ingroup prob_dists
 * The cumulative density function for a gamma distribution for y
 * with the specified shape and inverse scale parameters.
 *
 * @tparam T_y type of scalar
 * @tparam T_shape type of shape
 * @tparam T_inv_scale type of inverse scale
 * @param y A scalar variable.
 * @param alpha Shape parameter.
 * @param beta Inverse scale parameter.
 * @throw std::domain_error if y is not greater than or equal to 0.
 * @throw std::domain_error if alpha or beta is not greater than 0.
 */
template <typename T_y, typename T_shape, typename T_inv_scale>
return_type_t<T_y, T_shape, T_inv_scale> gamma_cdf(const T_y& y,
                                                   const T_shape& alpha,
                                                   const T_inv_scale& beta) {
  using T_partials_return = partials_return_t<T_y, T_shape, T_inv_scale>;
  using std::exp;
  static const char* function = "gamma_cdf";
  check_positive_finite(function, "Shape parameter", alpha);
  check_positive_finite(function, "Inverse scale parameter", beta);
  check_not_nan(function, "Random variable", y);
  check_nonnegative(function, "Random variable", y);
  check_consistent_sizes(function, "Random variable", y, "Shape parameter",
                         alpha, "Inverse scale parameter", beta);

  if (size_zero(y, alpha, beta)) {
    return 1.0;
  }

  T_partials_return P(1.0);
  operands_and_partials<T_y, T_shape, T_inv_scale> ops_partials(y, alpha, beta);

  scalar_seq_view<T_y> y_vec(y);
  scalar_seq_view<T_shape> alpha_vec(alpha);
  scalar_seq_view<T_inv_scale> beta_vec(beta);
  size_t size_y = stan::math::size(y);
  size_t size_alpha = stan::math::size(alpha);
  size_t N = max_size(y, alpha, beta);

  // Explicit return for extreme values
  // The gradients are technically ill-defined, but treated as zero
  for (size_t i = 0; i < size_y; i++) {
    if (value_of(y_vec[i]) == 0) {
      return ops_partials.build(0.0);
    }
  }

  VectorBuilder<!is_constant_all<T_y, T_shape, T_inv_scale>::value,
                T_partials_return, T_shape>
      tgamma_vec(size_alpha);
  VectorBuilder<!is_constant_all<T_shape>::value, T_partials_return, T_shape>
      digamma_vec(size_alpha);
  if (!is_constant_all<T_y, T_shape, T_inv_scale>::value) {
    for (size_t i = 0; i < size_alpha; i++) {
      const T_partials_return alpha_dbl = value_of(alpha_vec[i]);
      tgamma_vec[i] = tgamma(alpha_dbl);
      if (!is_constant_all<T_shape>::value) {
        digamma_vec[i] = digamma(alpha_dbl);
      }
    }
  }

  for (size_t n = 0; n < N; n++) {
    // Explicit results for extreme values
    // The gradients are technically ill-defined, but treated as zero
    if (value_of(y_vec[n]) == INFTY) {
      continue;
    }

    const T_partials_return y_dbl = value_of(y_vec[n]);
    const T_partials_return alpha_dbl = value_of(alpha_vec[n]);
    const T_partials_return beta_dbl = value_of(beta_vec[n]);
    const T_partials_return beta_y = beta_dbl * y_dbl;
    const T_partials_return Pn = gamma_p(alpha_dbl, beta_y);
    const T_partials_return rep_deriv = is_constant_all<T_y, T_inv_scale>::value
                                            ? 0
                                            : exp(-beta_y)
                                                  * pow(beta_y, alpha_dbl - 1)
                                                  / (tgamma_vec[n] * Pn);
    P *= Pn;

    if (!is_constant_all<T_y>::value) {
      ops_partials.edge1_.partials_[n] += beta_dbl * rep_deriv;
    }
    if (!is_constant_all<T_shape>::value) {
      ops_partials.edge2_.partials_[n]
          -= grad_reg_inc_gamma(alpha_dbl, beta_y, tgamma_vec[n],
                                digamma_vec[n])
             / Pn;
    }
    if (!is_constant_all<T_inv_scale>::value) {
      ops_partials.edge3_.partials_[n] += y_dbl * rep_deriv;
    }
  }

  if (!is_constant_all<T_y>::value) {
    for (size_t n = 0; n < size_y; ++n) {
      ops_partials.edge1_.partials_[n] *= P;
    }
  }
  if (!is_constant_all<T_shape>::value) {
    for (size_t n = 0; n < size_alpha; ++n) {
      ops_partials.edge2_.partials_[n] *= P;
    }
  }
  if (!is_constant_all<T_inv_scale>::value) {
    for (size_t n = 0; n < stan::math::size(beta); ++n) {
      ops_partials.edge3_.partials_[n] *= P;
    }
  }
  return ops_partials.build(P);
}

}  // namespace math
}  // namespace stan
#endif
