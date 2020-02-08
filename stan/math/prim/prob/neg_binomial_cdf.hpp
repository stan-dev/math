#ifndef STAN_MATH_PRIM_PROB_NEG_BINOMIAL_CDF_HPP
#define STAN_MATH_PRIM_PROB_NEG_BINOMIAL_CDF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/digamma.hpp>
#include <stan/math/prim/fun/inc_beta.hpp>
#include <stan/math/prim/fun/inc_beta_dda.hpp>
#include <stan/math/prim/fun/inc_beta_ddz.hpp>
#include <stan/math/prim/fun/inv.hpp>
#include <stan/math/prim/fun/max_size.hpp>
#include <stan/math/prim/fun/size.hpp>
#include <stan/math/prim/fun/size_zero.hpp>
#include <stan/math/prim/fun/square.hpp>
#include <stan/math/prim/fun/value_of.hpp>
#include <cmath>
#include <limits>

namespace stan {
namespace math {

template <typename T_n, typename T_shape, typename T_inv_scale>
return_type_t<T_shape, T_inv_scale> neg_binomial_cdf(const T_n& n,
                                                     const T_shape& alpha,
                                                     const T_inv_scale& beta) {
  static const char* function = "neg_binomial_cdf";
  using T_partials_return = partials_return_t<T_n, T_shape, T_inv_scale>;

  if (size_zero(n, alpha, beta)) {
    return 1.0;
  }

  T_partials_return P(1.0);

  check_positive_finite(function, "Shape parameter", alpha);
  check_positive_finite(function, "Inverse scale parameter", beta);
  check_consistent_sizes(function, "Failures variable", n, "Shape parameter",
                         alpha, "Inverse scale parameter", beta);

  scalar_seq_view<T_n> n_vec(n);
  scalar_seq_view<T_shape> alpha_vec(alpha);
  scalar_seq_view<T_inv_scale> beta_vec(beta);
  size_t size_alpha = size(alpha);
  size_t size_n_alpha = max_size(n, alpha);
  size_t max_size_seq_view = max_size(n, alpha, beta);

  operands_and_partials<T_shape, T_inv_scale> ops_partials(alpha, beta);

  // Explicit return for extreme values
  // The gradients are technically ill-defined, but treated as zero
  for (size_t i = 0; i < size(n); i++) {
    if (value_of(n_vec[i]) < 0) {
      return ops_partials.build(0.0);
    }
  }

  VectorBuilder<!is_constant_all<T_shape>::value, T_partials_return, T_shape>
      digamma_alpha_vec(size_alpha);
  VectorBuilder<!is_constant_all<T_shape>::value, T_partials_return, T_n,
                T_shape>
      digamma_sum_vec(size_n_alpha);

  if (!is_constant_all<T_shape>::value) {
    for (size_t i = 0; i < size_alpha; i++) {
      digamma_alpha_vec[i] = digamma(value_of(alpha_vec[i]));
    }
    for (size_t i = 0; i < size_n_alpha; i++) {
      const T_partials_return n_dbl = value_of(n_vec[i]);
      const T_partials_return alpha_dbl = value_of(alpha_vec[i]);
      digamma_sum_vec[i] = digamma(n_dbl + alpha_dbl + 1);
    }
  }

  for (size_t i = 0; i < max_size_seq_view; i++) {
    // Explicit results for extreme values
    // The gradients are technically ill-defined, but treated as zero
    if (value_of(n_vec[i]) == std::numeric_limits<int>::max()) {
      return ops_partials.build(1.0);
    }

    const T_partials_return n_dbl = value_of(n_vec[i]);
    const T_partials_return alpha_dbl = value_of(alpha_vec[i]);
    const T_partials_return beta_dbl = value_of(beta_vec[i]);
    const T_partials_return inv_beta_p1 = inv(beta_dbl + 1);
    const T_partials_return p_dbl = beta_dbl * inv_beta_p1;
    const T_partials_return d_dbl = square(inv_beta_p1);

    const T_partials_return P_i = inc_beta(alpha_dbl, n_dbl + 1.0, p_dbl);

    P *= P_i;

    if (!is_constant_all<T_shape>::value) {
      ops_partials.edge1_.partials_[i]
          += inc_beta_dda(alpha_dbl, n_dbl + 1, p_dbl, digamma_alpha_vec[i],
                          digamma_sum_vec[i])
             / P_i;
    }

    if (!is_constant_all<T_inv_scale>::value) {
      ops_partials.edge2_.partials_[i]
          += inc_beta_ddz(alpha_dbl, n_dbl + 1.0, p_dbl) * d_dbl / P_i;
    }
  }

  if (!is_constant_all<T_shape>::value) {
    for (size_t i = 0; i < size_alpha; ++i) {
      ops_partials.edge1_.partials_[i] *= P;
    }
  }

  if (!is_constant_all<T_inv_scale>::value) {
    for (size_t i = 0; i < size(beta); ++i) {
      ops_partials.edge2_.partials_[i] *= P;
    }
  }

  return ops_partials.build(P);
}

}  // namespace math
}  // namespace stan
#endif
