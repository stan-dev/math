#ifndef STAN_MATH_PRIM_PROB_NEG_BINOMIAL_LCDF_HPP
#define STAN_MATH_PRIM_PROB_NEG_BINOMIAL_LCDF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/size_zero.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/value_of.hpp>
#include <stan/math/prim/fun/digamma.hpp>
#include <stan/math/prim/fun/beta.hpp>
#include <stan/math/prim/fun/grad_reg_inc_beta.hpp>
#include <stan/math/prim/fun/inc_beta.hpp>
#include <cmath>
#include <limits>

namespace stan {
namespace math {

template <typename T_n, typename T_shape, typename T_inv_scale>
return_type_t<T_shape, T_inv_scale> neg_binomial_lcdf(const T_n& n,
                                                      const T_shape& alpha,
                                                      const T_inv_scale& beta) {
  static const char* function = "neg_binomial_lcdf";
  using T_partials_return = partials_return_t<T_n, T_shape, T_inv_scale>;

  if (size_zero(n, alpha, beta)) {
    return 0.0;
  }

  T_partials_return P(0.0);

  check_positive_finite(function, "Shape parameter", alpha);
  check_positive_finite(function, "Inverse scale parameter", beta);
  check_consistent_sizes(function, "Failures variable", n, "Shape parameter",
                         alpha, "Inverse scale parameter", beta);

  scalar_seq_view<T_n> n_vec(n);
  scalar_seq_view<T_shape> alpha_vec(alpha);
  scalar_seq_view<T_inv_scale> beta_vec(beta);
  size_t max_size_seq_view = max_size(n, alpha, beta);

  using std::exp;
  using std::log;
  using std::pow;

  operands_and_partials<T_shape, T_inv_scale> ops_partials(alpha, beta);

  // Explicit return for extreme values
  // The gradients are technically ill-defined, but treated as zero
  for (size_t i = 0; i < size(n); i++) {
    if (value_of(n_vec[i]) < 0) {
      return ops_partials.build(negative_infinity());
    }
  }

  VectorBuilder<!is_constant_all<T_shape>::value, T_partials_return, T_shape>
      digammaN_vec(size(alpha));
  VectorBuilder<!is_constant_all<T_shape>::value, T_partials_return, T_shape>
      digammaAlpha_vec(size(alpha));
  VectorBuilder<!is_constant_all<T_shape>::value, T_partials_return, T_shape>
      digammaSum_vec(size(alpha));

  if (!is_constant_all<T_shape>::value) {
    for (size_t i = 0; i < size(alpha); i++) {
      const T_partials_return n_dbl = value_of(n_vec[i]);
      const T_partials_return alpha_dbl = value_of(alpha_vec[i]);

      digammaN_vec[i] = digamma(n_dbl + 1);
      digammaAlpha_vec[i] = digamma(alpha_dbl);
      digammaSum_vec[i] = digamma(n_dbl + alpha_dbl + 1);
    }
  }

  for (size_t i = 0; i < max_size_seq_view; i++) {
    // Explicit results for extreme values
    // The gradients are technically ill-defined, but treated as zero
    if (value_of(n_vec[i]) == std::numeric_limits<int>::max()) {
      return ops_partials.build(0.0);
    }

    const T_partials_return n_dbl = value_of(n_vec[i]);
    const T_partials_return alpha_dbl = value_of(alpha_vec[i]);
    const T_partials_return beta_dbl = value_of(beta_vec[i]);
    const T_partials_return p_dbl = beta_dbl / (1.0 + beta_dbl);
    const T_partials_return d_dbl = 1.0 / ((1.0 + beta_dbl) * (1.0 + beta_dbl));
    const T_partials_return Pi = inc_beta(alpha_dbl, n_dbl + 1.0, p_dbl);
    const T_partials_return beta_func = stan::math::beta(n_dbl + 1, alpha_dbl);

    P += log(Pi);

    if (!is_constant_all<T_shape>::value) {
      T_partials_return g1 = 0;
      T_partials_return g2 = 0;

      grad_reg_inc_beta(g1, g2, alpha_dbl, n_dbl + 1, p_dbl,
                        digammaAlpha_vec[i], digammaN_vec[i], digammaSum_vec[i],
                        beta_func);
      ops_partials.edge1_.partials_[i] += g1 / Pi;
    }
    if (!is_constant_all<T_inv_scale>::value) {
      ops_partials.edge2_.partials_[i] += d_dbl * pow(1 - p_dbl, n_dbl)
                                          * pow(p_dbl, alpha_dbl - 1)
                                          / beta_func / Pi;
    }
  }

  return ops_partials.build(P);
}

}  // namespace math
}  // namespace stan
#endif
