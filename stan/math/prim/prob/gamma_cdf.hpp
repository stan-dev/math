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
#include <stan/math/prim/fun/scalar_seq_view.hpp>
#include <stan/math/prim/fun/size.hpp>
#include <stan/math/prim/fun/size_zero.hpp>
#include <stan/math/prim/fun/tgamma.hpp>
#include <stan/math/prim/fun/value_of.hpp>
#include <boost/random/gamma_distribution.hpp>
#include <boost/random/variate_generator.hpp>
#include <stan/math/prim/functor/partials_propagator.hpp>
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
 * @throw std::domain_error if alpha is not greater than 0.
 * @throw std::domain_error if beta is not greater than 0.
 * @throw std::domain_error if y is not greater than or equal to 0.
 */
template <typename T_y, typename T_shape, typename T_inv_scale>
return_type_t<T_y, T_shape, T_inv_scale> gamma_cdf(const T_y& y,
                                                   const T_shape& alpha,
                                                   const T_inv_scale& beta) {
  using T_partials_return = partials_return_t<T_y, T_shape, T_inv_scale>;
  using T_y_ref = ref_type_t<T_y>;
  using T_alpha_ref = ref_type_t<T_shape>;
  using T_beta_ref = ref_type_t<T_inv_scale>;
  using std::exp;
  static const char* function = "gamma_cdf";
  check_consistent_sizes(function, "Random variable", y, "Shape parameter",
                         alpha, "Inverse scale parameter", beta);
  T_y_ref y_ref = y;
  T_alpha_ref alpha_ref = alpha;
  T_beta_ref beta_ref = beta;
  check_positive_finite(function, "Shape parameter", alpha_ref);
  check_positive_finite(function, "Inverse scale parameter", beta_ref);
  check_nonnegative(function, "Random variable", y_ref);

  if (size_zero(y, alpha, beta)) {
    return 1.0;
  }

  T_partials_return P(1.0);
  auto ops_partials = make_partials_propagator(y_ref, alpha_ref, beta_ref);

  scalar_seq_view<T_y_ref> y_vec(y_ref);
  scalar_seq_view<T_alpha_ref> alpha_vec(alpha_ref);
  scalar_seq_view<T_beta_ref> beta_vec(beta_ref);
  size_t N = max_size(y, alpha, beta);

  // Explicit return for extreme values
  // The gradients are technically ill-defined, but treated as zero
  for (size_t i = 0; i < stan::math::size(y); i++) {
    if (y_vec.val(i) == 0) {
      return ops_partials.build(0.0);
    }
  }

  using std::exp;
  using std::pow;

  VectorBuilder<!is_constant_all<T_shape>::value, T_partials_return, T_shape>
      gamma_vec(math::size(alpha));
  VectorBuilder<!is_constant_all<T_shape>::value, T_partials_return, T_shape>
      digamma_vec(math::size(alpha));

  if (!is_constant_all<T_shape>::value) {
    for (size_t i = 0; i < stan::math::size(alpha); i++) {
      const T_partials_return alpha_dbl = alpha_vec.val(i);
      gamma_vec[i] = tgamma(alpha_dbl);
      digamma_vec[i] = digamma(alpha_dbl);
    }
  }

  for (size_t n = 0; n < N; n++) {
    // Explicit results for extreme values
    // The gradients are technically ill-defined, but treated as zero
    if (y_vec.val(n) == INFTY) {
      continue;
    }

    const T_partials_return y_dbl = y_vec.val(n);
    const T_partials_return alpha_dbl = alpha_vec.val(n);
    const T_partials_return beta_dbl = beta_vec.val(n);

    const T_partials_return Pn = gamma_p(alpha_dbl, beta_dbl * y_dbl);

    P *= Pn;

    if (!is_constant_all<T_y>::value) {
      partials<0>(ops_partials)[n] += beta_dbl * exp(-beta_dbl * y_dbl)
                                      * pow(beta_dbl * y_dbl, alpha_dbl - 1)
                                      / tgamma(alpha_dbl) / Pn;
    }
    if (!is_constant_all<T_shape>::value) {
      partials<1>(ops_partials)[n]
          -= grad_reg_inc_gamma(alpha_dbl, beta_dbl * y_dbl, gamma_vec[n],
                                digamma_vec[n])
             / Pn;
    }
    if (!is_constant_all<T_inv_scale>::value) {
      partials<2>(ops_partials)[n] += y_dbl * exp(-beta_dbl * y_dbl)
                                      * pow(beta_dbl * y_dbl, alpha_dbl - 1)
                                      / tgamma(alpha_dbl) / Pn;
    }
  }

  if (!is_constant_all<T_y>::value) {
    for (size_t n = 0; n < stan::math::size(y); ++n) {
      partials<0>(ops_partials)[n] *= P;
    }
  }
  if (!is_constant_all<T_shape>::value) {
    for (size_t n = 0; n < stan::math::size(alpha); ++n) {
      partials<1>(ops_partials)[n] *= P;
    }
  }
  if (!is_constant_all<T_inv_scale>::value) {
    for (size_t n = 0; n < stan::math::size(beta); ++n) {
      partials<2>(ops_partials)[n] *= P;
    }
  }
  return ops_partials.build(P);
}

}  // namespace math
}  // namespace stan
#endif
