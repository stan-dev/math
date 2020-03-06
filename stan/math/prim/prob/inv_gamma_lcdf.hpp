#ifndef STAN_MATH_PRIM_PROB_INV_GAMMA_LCDF_HPP
#define STAN_MATH_PRIM_PROB_INV_GAMMA_LCDF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/digamma.hpp>
#include <stan/math/prim/fun/exp.hpp>
#include <stan/math/prim/fun/gamma_q.hpp>
#include <stan/math/prim/fun/grad_reg_inc_gamma.hpp>
#include <stan/math/prim/fun/inv.hpp>
#include <stan/math/prim/fun/log.hpp>
#include <stan/math/prim/fun/max_size.hpp>
#include <stan/math/prim/fun/size.hpp>
#include <stan/math/prim/fun/size_zero.hpp>
#include <stan/math/prim/fun/tgamma.hpp>
#include <stan/math/prim/fun/value_of.hpp>
#include <cmath>

namespace stan {
namespace math {

template <typename T_y, typename T_shape, typename T_scale>
return_type_t<T_y, T_shape, T_scale> inv_gamma_lcdf(const T_y& y,
                                                    const T_shape& alpha,
                                                    const T_scale& beta) {
  using T_partials_return = partials_return_t<T_y, T_shape, T_scale>;
  using std::exp;
  using std::log;
  using std::pow;
  static const char* function = "inv_gamma_lcdf";
  check_positive_finite(function, "Shape parameter", alpha);
  check_positive_finite(function, "Scale parameter", beta);
  check_not_nan(function, "Random variable", y);
  check_nonnegative(function, "Random variable", y);
  check_consistent_sizes(function, "Random variable", y, "Shape parameter",
                         alpha, "Scale Parameter", beta);

  if (size_zero(y, alpha, beta)) {
    return 0;
  }

  T_partials_return P(0.0);
  operands_and_partials<T_y, T_shape, T_scale> ops_partials(y, alpha, beta);

  scalar_seq_view<T_y> y_vec(y);
  scalar_seq_view<T_shape> alpha_vec(alpha);
  scalar_seq_view<T_scale> beta_vec(beta);
  size_t size_y = stan::math::size(y);
  size_t size_alpha = stan::math::size(alpha);
  size_t N = max_size(y, alpha, beta);

  // Explicit return for extreme values
  // The gradients are technically ill-defined, but treated as zero
  for (size_t i = 0; i < size_y; i++) {
    if (value_of(y_vec[i]) == 0) {
      return ops_partials.build(NEGATIVE_INFTY);
    }
  }

  VectorBuilder<true, T_partials_return, T_y> inv_y(size_y);
  for (size_t i = 0; i < size_y; i++) {
    inv_y[i] = inv(value_of(y_vec[i]));
  }

  VectorBuilder<!is_constant_all<T_y, T_shape, T_scale>::value,
                T_partials_return, T_shape>
      tgamma_vec(size_alpha);
  VectorBuilder<!is_constant_all<T_shape>::value, T_partials_return, T_shape>
      digamma_vec(size_alpha);
  if (!is_constant_all<T_y, T_shape, T_scale>::value) {
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

    const T_partials_return alpha_dbl = value_of(alpha_vec[n]);
    const T_partials_return beta_over_y = value_of(beta_vec[n]) * inv_y[n];
    const T_partials_return Pn = gamma_q(alpha_dbl, beta_over_y);
    const T_partials_return rep_deriv
        = is_constant_all<T_y, T_scale>::value
              ? 0
              : inv_y[n] * exp(-beta_over_y) * pow(beta_over_y, alpha_dbl - 1)
                    / (tgamma_vec[n] * Pn);

    P += log(Pn);

    if (!is_constant_all<T_y>::value) {
      ops_partials.edge1_.partials_[n] += beta_over_y * rep_deriv;
    }
    if (!is_constant_all<T_shape>::value) {
      ops_partials.edge2_.partials_[n]
          += grad_reg_inc_gamma(alpha_dbl, beta_over_y, tgamma_vec[n],
                                digamma_vec[n])
             / Pn;
    }
    if (!is_constant_all<T_scale>::value) {
      ops_partials.edge3_.partials_[n] -= rep_deriv;
    }
  }
  return ops_partials.build(P);
}

}  // namespace math
}  // namespace stan
#endif
