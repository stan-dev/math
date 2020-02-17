#ifndef STAN_MATH_PRIM_PROB_INV_GAMMA_CDF_HPP
#define STAN_MATH_PRIM_PROB_INV_GAMMA_CDF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/digamma.hpp>
#include <stan/math/prim/fun/exp.hpp>
#include <stan/math/prim/fun/gamma_q.hpp>
#include <stan/math/prim/fun/grad_reg_inc_gamma.hpp>
#include <stan/math/prim/fun/max_size.hpp>
#include <stan/math/prim/fun/size.hpp>
#include <stan/math/prim/fun/size_zero.hpp>
#include <stan/math/prim/fun/tgamma.hpp>
#include <stan/math/prim/fun/value_of.hpp>
#include <cmath>

namespace stan {
namespace math {

/** \ingroup prob_dists
 * The CDF of an inverse gamma density for y with the specified
 * shape and scale parameters. y, shape, and scale parameters must
 * be greater than 0.
 *
 * @param y A scalar variable.
 * @param alpha Shape parameter.
 * @param beta Scale parameter.
 * @throw std::domain_error if alpha is not greater than 0.
 * @throw std::domain_error if beta is not greater than 0.
 * @throw std::domain_error if y is not greater than 0.
 * @tparam T_y Type of scalar.
 * @tparam T_shape Type of shape.
 * @tparam T_scale Type of scale.
 */

template <typename T_y, typename T_shape, typename T_scale>
return_type_t<T_y, T_shape, T_scale> inv_gamma_cdf(const T_y& y,
                                                   const T_shape& alpha,
                                                   const T_scale& beta) {
  using T_partials_return = partials_return_t<T_y, T_shape, T_scale>;

  if (size_zero(y, alpha, beta)) {
    return 1.0;
  }

  static const char* function = "inv_gamma_cdf";

  T_partials_return P(1.0);

  check_positive_finite(function, "Shape parameter", alpha);
  check_positive_finite(function, "Scale parameter", beta);
  check_not_nan(function, "Random variable", y);
  check_nonnegative(function, "Random variable", y);
  check_consistent_sizes(function, "Random variable", y, "Shape parameter",
                         alpha, "Scale Parameter", beta);

  scalar_seq_view<T_y> y_vec(y);
  scalar_seq_view<T_shape> alpha_vec(alpha);
  scalar_seq_view<T_scale> beta_vec(beta);
  size_t N = max_size(y, alpha, beta);

  operands_and_partials<T_y, T_shape, T_scale> ops_partials(y, alpha, beta);

  // Explicit return for extreme values
  // The gradients are technically ill-defined, but treated as zero
  for (size_t i = 0; i < size(y); i++) {
    if (value_of(y_vec[i]) == 0) {
      return ops_partials.build(0.0);
    }
  }

  using std::exp;
  using std::pow;

  VectorBuilder<!is_constant_all<T_shape>::value, T_partials_return, T_shape>
      gamma_vec(size(alpha));
  VectorBuilder<!is_constant_all<T_shape>::value, T_partials_return, T_shape>
      digamma_vec(size(alpha));

  if (!is_constant_all<T_shape>::value) {
    for (size_t i = 0; i < size(alpha); i++) {
      const T_partials_return alpha_dbl = value_of(alpha_vec[i]);
      gamma_vec[i] = tgamma(alpha_dbl);
      digamma_vec[i] = digamma(alpha_dbl);
    }
  }

  for (size_t n = 0; n < N; n++) {
    // Explicit results for extreme values
    // The gradients are technically ill-defined, but treated as zero
    if (value_of(y_vec[n]) == INFTY) {
      continue;
    }

    const T_partials_return y_dbl = value_of(y_vec[n]);
    const T_partials_return y_inv_dbl = 1.0 / y_dbl;
    const T_partials_return alpha_dbl = value_of(alpha_vec[n]);
    const T_partials_return beta_dbl = value_of(beta_vec[n]);

    const T_partials_return Pn = gamma_q(alpha_dbl, beta_dbl * y_inv_dbl);

    P *= Pn;

    if (!is_constant_all<T_y>::value) {
      ops_partials.edge1_.partials_[n]
          += beta_dbl * y_inv_dbl * y_inv_dbl * exp(-beta_dbl * y_inv_dbl)
             * pow(beta_dbl * y_inv_dbl, alpha_dbl - 1) / tgamma(alpha_dbl)
             / Pn;
    }
    if (!is_constant_all<T_shape>::value) {
      ops_partials.edge2_.partials_[n]
          += grad_reg_inc_gamma(alpha_dbl, beta_dbl * y_inv_dbl, gamma_vec[n],
                                digamma_vec[n])
             / Pn;
    }
    if (!is_constant_all<T_scale>::value) {
      ops_partials.edge3_.partials_[n]
          += -y_inv_dbl * exp(-beta_dbl * y_inv_dbl)
             * pow(beta_dbl * y_inv_dbl, alpha_dbl - 1) / tgamma(alpha_dbl)
             / Pn;
    }
  }

  if (!is_constant_all<T_y>::value) {
    for (size_t n = 0; n < size(y); ++n) {
      ops_partials.edge1_.partials_[n] *= P;
    }
  }
  if (!is_constant_all<T_shape>::value) {
    for (size_t n = 0; n < size(alpha); ++n) {
      ops_partials.edge2_.partials_[n] *= P;
    }
  }
  if (!is_constant_all<T_scale>::value) {
    for (size_t n = 0; n < size(beta); ++n) {
      ops_partials.edge3_.partials_[n] *= P;
    }
  }
  return ops_partials.build(P);
}

}  // namespace math
}  // namespace stan
#endif
