#ifndef STAN_MATH_PRIM_PROB_GAMMA_LCCDF_HPP
#define STAN_MATH_PRIM_PROB_GAMMA_LCCDF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/digamma.hpp>
#include <stan/math/prim/fun/gamma_q.hpp>
#include <stan/math/prim/fun/grad_reg_inc_gamma.hpp>
#include <stan/math/prim/fun/size_zero.hpp>
#include <stan/math/prim/fun/tgamma.hpp>
#include <stan/math/prim/fun/value_of.hpp>
#include <cmath>

namespace stan {
namespace math {

template <typename T_y, typename T_shape, typename T_inv_scale>
return_type_t<T_y, T_shape, T_inv_scale> gamma_lccdf(const T_y& y,
                                                     const T_shape& alpha,
                                                     const T_inv_scale& beta) {
  if (size_zero(y, alpha, beta)) {
    return 0.0;
  }

  using T_partials_return = partials_return_t<T_y, T_shape, T_inv_scale>;

  static const char* function = "gamma_lccdf";

  T_partials_return P(0.0);

  check_positive_finite(function, "Shape parameter", alpha);
  check_positive_finite(function, "Inverse scale parameter", beta);
  check_not_nan(function, "Random variable", y);
  check_nonnegative(function, "Random variable", y);
  check_consistent_sizes(function, "Random variable", y, "Shape parameter",
                         alpha, "Inverse scale parameter", beta);

  scalar_seq_view<T_y> y_vec(y);
  scalar_seq_view<T_shape> alpha_vec(alpha);
  scalar_seq_view<T_inv_scale> beta_vec(beta);
  size_t N = max_size(y, alpha, beta);

  operands_and_partials<T_y, T_shape, T_inv_scale> ops_partials(y, alpha, beta);

  // Explicit return for extreme values
  // The gradients are technically ill-defined, but treated as zero
  for (size_t i = 0; i < size(y); i++) {
    if (value_of(y_vec[i]) == 0) {
      return ops_partials.build(0.0);
    }
  }

  using std::exp;
  using std::log;
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
      return ops_partials.build(negative_infinity());
    }

    const T_partials_return y_dbl = value_of(y_vec[n]);
    const T_partials_return alpha_dbl = value_of(alpha_vec[n]);
    const T_partials_return beta_dbl = value_of(beta_vec[n]);

    const T_partials_return Pn = gamma_q(alpha_dbl, beta_dbl * y_dbl);

    P += log(Pn);

    if (!is_constant_all<T_y>::value) {
      ops_partials.edge1_.partials_[n] -= beta_dbl * exp(-beta_dbl * y_dbl)
                                          * pow(beta_dbl * y_dbl, alpha_dbl - 1)
                                          / tgamma(alpha_dbl) / Pn;
    }
    if (!is_constant_all<T_shape>::value) {
      ops_partials.edge2_.partials_[n]
          += grad_reg_inc_gamma(alpha_dbl, beta_dbl * y_dbl, gamma_vec[n],
                                digamma_vec[n])
             / Pn;
    }
    if (!is_constant_all<T_inv_scale>::value) {
      ops_partials.edge3_.partials_[n] -= y_dbl * exp(-beta_dbl * y_dbl)
                                          * pow(beta_dbl * y_dbl, alpha_dbl - 1)
                                          / tgamma(alpha_dbl) / Pn;
    }
  }
  return ops_partials.build(P);
}

}  // namespace math
}  // namespace stan

#endif
